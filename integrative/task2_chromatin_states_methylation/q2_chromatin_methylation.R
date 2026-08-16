library(GenomicRanges)
library(dplyr)
library(ggplot2)

kidney_states <- read.table("/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/segmentation/kidney_15_dense.bed", skip = 1, sep = "\t")[,1:4]
liver_states <- read.table("/vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/segmentation/liver_15_dense.bed", skip = 1, sep = "\t")[,1:4]
colnames(kidney_states) <- c("chr", "start", "end", "state")
colnames(liver_states) <- c("chr", "start", "end", "state")

kidney_meth <- read.table("/vol/COMPEPIWS/groups/shared/WGBS/wgbs1/signal/rnbeads_kidney.bedGraph", sep = "\t", skip = 1)
colnames(kidney_meth) <- c("chr", "start", "end", "meth")
liver_meth <- read.table("/vol/COMPEPIWS/groups/shared/WGBS/wgbs1/signal/rnbeads_liver.bedGraph", sep = "\t", skip = 1)
colnames(liver_meth) <- c("chr", "start", "end", "meth")

kidney_states_gr <- GRanges(kidney_states$chr, IRanges(kidney_states$start + 1, kidney_states$end), state = kidney_states$state)
liver_states_gr <- GRanges(liver_states$chr, IRanges(liver_states$start + 1, liver_states$end), state = liver_states$state)
kidney_meth_gr <- GRanges(kidney_meth$chr, IRanges(kidney_meth$start + 1, kidney_meth$end), meth = kidney_meth$meth)
liver_meth_gr <- GRanges(liver_meth$chr, IRanges(liver_meth$start + 1, liver_meth$end), meth = liver_meth$meth)

hits_kidney <- findOverlaps(kidney_meth_gr, kidney_states_gr)
kidney_df <- data.frame(meth = kidney_meth_gr$meth[queryHits(hits_kidney)], state = kidney_states_gr$state[subjectHits(hits_kidney)], sample = "kidney")

hits_liver <- findOverlaps(liver_meth_gr, liver_states_gr)
liver_df <- data.frame(meth = liver_meth_gr$meth[queryHits(hits_liver)], state = liver_states_gr$state[subjectHits(hits_liver)], sample = "liver")

combined_df <- rbind(kidney_df, liver_df)

avg_meth <- combined_df %>% group_by(sample, state) %>% summarise(mean_meth = mean(meth, na.rm = TRUE))
print(avg_meth, n = 40)
write.csv(avg_meth, "avg_methylation_per_state.csv", row.names = FALSE)

png("methylation_by_state_boxplot.png", width = 1200, height = 900)
ggplot(combined_df, aes(x = state, y = meth, fill = sample)) + geom_boxplot(outlier.size = 0.3) + facet_wrap(~sample, ncol = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
