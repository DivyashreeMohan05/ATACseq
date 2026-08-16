library(GenomicRanges)
library(dplyr)
library(ComplexHeatmap)

dmr <- read.table("dmr_closest_gene.bed", sep = "\t")
colnames(dmr) <- c("chr","start","end","score","g_chr","g_start","g_end","transcript_id","x1","strand","x2","x3","x4","x5","x6","distance")
dmr <- dmr[, c("chr","start","end","score","transcript_id","distance")]

map <- read.table("transcript_to_gene.tsv", sep = "\t")
colnames(map) <- c("transcript_id","gene")
dmr <- merge(dmr, map, by = "transcript_id")

kidney_meth <- read.table("/vol/COMPEPIWS/groups/shared/WGBS/wgbs1/signal/rnbeads_kidney.bedGraph", sep = "\t", skip = 1)
colnames(kidney_meth) <- c("chr","start","end","meth")
liver_meth <- read.table("/vol/COMPEPIWS/groups/shared/WGBS/wgbs1/signal/rnbeads_liver.bedGraph", sep = "\t", skip = 1)
colnames(liver_meth) <- c("chr","start","end","meth")

dmr_gr <- GRanges(dmr$chr, IRanges(dmr$start + 1, dmr$end))
kidney_gr <- GRanges(kidney_meth$chr, IRanges(kidney_meth$start + 1, kidney_meth$end), meth = kidney_meth$meth)
liver_gr <- GRanges(liver_meth$chr, IRanges(liver_meth$start + 1, liver_meth$end), meth = liver_meth$meth)

hk <- findOverlaps(dmr_gr, kidney_gr)
hl <- findOverlaps(dmr_gr, liver_gr)

dmr$kidney_meth <- NA
dmr$liver_meth <- NA
dmr$kidney_meth[unique(queryHits(hk))] <- sapply(unique(queryHits(hk)), function(i) mean(kidney_gr$meth[subjectHits(hk)[queryHits(hk)==i]]))
dmr$liver_meth[unique(queryHits(hl))] <- sapply(unique(queryHits(hl)), function(i) mean(liver_gr$meth[subjectHits(hl)[queryHits(hl)==i]]))

cpm <- read.csv("/vol/COMPEPIWS/groups/shared/RNA-seq/rnaseq1/DEGs/CPM_expression.csv")
cpm$kidney_cpm <- rowMeans(cpm[, c("kidney_14_5_RNA_1","kidney_14_5_RNA_2","kidney_15_5_RNA_1","kidney_15_5_RNA_2")])
cpm$liver_cpm <- rowMeans(cpm[, c("liver_14_5_RNA_1","liver_14_5_RNA_2","liver_15_5_RNA_1","liver_15_5_RNA_2")])
cpm_small <- cpm[, c("gene_ID","kidney_cpm","liver_cpm")]
colnames(cpm_small) <- c("gene","kidney_cpm","liver_cpm")

gene_summary <- dmr %>%
  group_by(gene) %>%
  summarise(mean_kidney_meth = mean(kidney_meth, na.rm = TRUE),
            mean_liver_meth = mean(liver_meth, na.rm = TRUE),
            n_dmrs = n())
gene_summary <- merge(gene_summary, cpm_small, by = "gene")

write.csv(gene_summary, "gene_summary.csv", row.names = FALSE)

top_genes <- gene_summary %>% arrange(desc(n_dmrs)) %>% head(30)
meth_matrix <- as.matrix(top_genes[, c("mean_kidney_meth","mean_liver_meth")])
rownames(meth_matrix) <- top_genes$gene
colnames(meth_matrix) <- c("Kidney","Liver")

row_order_methylation <- hclust(dist(meth_matrix))$order
gene_order <- rownames(meth_matrix)[row_order_methylation]

meth_matrix_ordered <- meth_matrix[gene_order, ]

cpm_matrix <- as.matrix(log2(top_genes[, c("kidney_cpm","liver_cpm")] + 1))
rownames(cpm_matrix) <- top_genes$gene
colnames(cpm_matrix) <- c("Kidney","Liver")
cpm_matrix_ordered <- cpm_matrix[gene_order, ]

png("heatmap_methylation.png", width = 800, height = 1000)
Heatmap(meth_matrix_ordered, name = "Methylation", cluster_columns = FALSE, cluster_rows = FALSE,
        column_title = "Mean Methylation per Group (Top 30 DMR genes)")
dev.off()

png("heatmap_cpm.png", width = 800, height = 1000)
Heatmap(cpm_matrix_ordered, name = "log2(CPM+1)", cluster_columns = FALSE, cluster_rows = FALSE,
        column_title = "Mean log2(CPM+1) per Group (Top 30 DMR genes)")
dev.off()

cor.test(gene_summary$mean_kidney_meth, gene_summary$kidney_cpm, method = "spearman")
cor.test(gene_summary$mean_liver_meth, gene_summary$liver_cpm, method = "spearman")
