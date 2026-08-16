# Integrative Task 4
# Kidney-liver DMR overlap with ChromHMM Het states
# and differential ATAC-seq peaks

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

dmr_file <- paste0(
  "/vol/COMPEPIWS/groups/shared/WGBS/",
  "wgbs1/differential/diff_meth_table.csv"
)

kidney_chromhmm_file <- paste0(
  "/vol/COMPEPIWS/groups/shared/ChIP-seq/",
  "chipseq1/segmentation/kidney_15_segments.bed"
)

liver_chromhmm_file <- paste0(
  "/vol/COMPEPIWS/groups/shared/ChIP-seq/",
  "chipseq1/segmentation/liver_15_segments.bed"
)

kidney_atac_file <- paste0(
  "/vol/COMPEPIWS/groups/shared/ATAC-seq/",
  "atacseq1/differential/kidney_specific_peaks.bed"
)

liver_atac_file <- paste0(
  "/vol/COMPEPIWS/groups/shared/ATAC-seq/",
  "atacseq1/differential/liver_specific_peaks.bed"
)

# ------------------------------------------------------------
# 2. Load and filter kidney-liver DMRs
# ------------------------------------------------------------

dmr_all <- read.csv(
  dmr_file,
  stringsAsFactors = FALSE,
  check.names = TRUE
)

# Define significant DMRs using FDR < 0.05
dmr_sig <- dmr_all[
  !is.na(dmr_all$comb.p.adj.fdr) &
    dmr_all$comb.p.adj.fdr < 0.05,
]

# WGBS CSV coordinates are treated as 1-based closed coordinates
dmr_gr <- GRanges(
  seqnames = dmr_sig$Chromosome,
  ranges = IRanges(
    start = dmr_sig$Start,
    end = dmr_sig$End
  )
)

mcols(dmr_gr)$DMR_ID <- if ("name" %in% colnames(dmr_sig)) {
  dmr_sig$name
} else {
  paste0("DMR_", seq_len(nrow(dmr_sig)))
}

mcols(dmr_gr)$meth_diff <- dmr_sig$mean.mean.diff
mcols(dmr_gr)$FDR <- dmr_sig$comb.p.adj.fdr

# ------------------------------------------------------------
# 3. Function to load BED files
# ------------------------------------------------------------

read_bed_granges <- function(file, has_state = FALSE) {

  bed <- read.delim(
    file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  if (ncol(bed) < 3) {
    stop("BED file has fewer than three columns: ", file)
  }

  if (has_state) {
    if (ncol(bed) < 4) {
      stop("ChromHMM file has no state column: ", file)
    }

    colnames(bed)[1:4] <- c(
      "chromosome",
      "start0",
      "end",
      "state"
    )
  } else {
    colnames(bed)[1:3] <- c(
      "chromosome",
      "start0",
      "end"
    )
  }

  # BED is 0-based half-open.
  # GRanges is 1-based closed.
  gr <- GRanges(
    seqnames = bed$chromosome,
    ranges = IRanges(
      start = bed$start0 + 1,
      end = bed$end
    )
  )

  if (has_state) {
    mcols(gr)$state <- bed$state
  }

  gr
}

# ------------------------------------------------------------
# 4. Load ChromHMM segmentations
# ------------------------------------------------------------

kidney_chromhmm_gr <- read_bed_granges(
  kidney_chromhmm_file,
  has_state = TRUE
)

liver_chromhmm_gr <- read_bed_granges(
  liver_chromhmm_file,
  has_state = TRUE
)

# Treat both ChromHMM Het subclasses as heterochromatin
kidney_het_gr <- kidney_chromhmm_gr[
  kidney_chromhmm_gr$state %in% c("Het_P", "Het_S")
]

liver_het_gr <- liver_chromhmm_gr[
  liver_chromhmm_gr$state %in% c("Het_P", "Het_S")
]

# ------------------------------------------------------------
# 5. Load differential ATAC peaks
# ------------------------------------------------------------

kidney_atac_gr <- read_bed_granges(
  kidney_atac_file,
  has_state = FALSE
)

liver_atac_gr <- read_bed_granges(
  liver_atac_file,
  has_state = FALSE
)

# ------------------------------------------------------------
# 6. Calculate overlaps
# ------------------------------------------------------------

dmr_in_kidney_het <- countOverlaps(
  dmr_gr,
  kidney_het_gr
) > 0

dmr_in_liver_het <- countOverlaps(
  dmr_gr,
  liver_het_gr
) > 0

dmr_in_any_het <- (
  dmr_in_kidney_het |
    dmr_in_liver_het
)

dmr_in_both_het <- (
  dmr_in_kidney_het &
    dmr_in_liver_het
)

dmr_overlap_kidney_atac <- countOverlaps(
  dmr_gr,
  kidney_atac_gr
) > 0

dmr_overlap_liver_atac <- countOverlaps(
  dmr_gr,
  liver_atac_gr
) > 0

dmr_overlap_any_atac <- (
  dmr_overlap_kidney_atac |
    dmr_overlap_liver_atac
)

dmr_overlap_both_atac <- (
  dmr_overlap_kidney_atac &
    dmr_overlap_liver_atac
)

# ------------------------------------------------------------
# 7. Prepare summary table
# ------------------------------------------------------------

total_dmrs <- length(dmr_gr)

summary_table <- data.frame(
  category = c(
    "Total significant DMRs",
    "DMRs overlapping kidney Het",
    "DMRs overlapping liver Het",
    "DMRs overlapping Het in either tissue",
    "DMRs overlapping Het in both tissues",
    "DMRs overlapping kidney-specific ATAC peaks",
    "DMRs overlapping liver-specific ATAC peaks",
    "DMRs overlapping any differential ATAC peak",
    "DMRs overlapping both ATAC peak sets"
  ),
  count = c(
    total_dmrs,
    sum(dmr_in_kidney_het),
    sum(dmr_in_liver_het),
    sum(dmr_in_any_het),
    sum(dmr_in_both_het),
    sum(dmr_overlap_kidney_atac),
    sum(dmr_overlap_liver_atac),
    sum(dmr_overlap_any_atac),
    sum(dmr_overlap_both_atac)
  )
)

summary_table$percentage_of_DMRs <- round(
  100 * summary_table$count / total_dmrs,
  2
)

# Total DMR row should be 100%
summary_table$percentage_of_DMRs[1] <- 100

# ------------------------------------------------------------
# 8. Prepare DMR-level result table
# ------------------------------------------------------------

dmr_results <- dmr_sig

dmr_results$overlap_kidney_Het <- dmr_in_kidney_het
dmr_results$overlap_liver_Het <- dmr_in_liver_het
dmr_results$overlap_any_Het <- dmr_in_any_het
dmr_results$overlap_both_Het <- dmr_in_both_het

dmr_results$overlap_kidney_specific_ATAC <- (
  dmr_overlap_kidney_atac
)

dmr_results$overlap_liver_specific_ATAC <- (
  dmr_overlap_liver_atac
)

dmr_results$overlap_any_differential_ATAC <- (
  dmr_overlap_any_atac
)

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

write.table(
  summary_table,
  file = "task4_overlap_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  dmr_results,
  file = "task4_DMR_overlap_details.csv",
  row.names = FALSE
)

# Save significant DMRs as proper BED3 coordinates
dmr_bed <- data.frame(
  chromosome = as.character(seqnames(dmr_gr)),
  start0 = start(dmr_gr) - 1,
  end = end(dmr_gr)
)

write.table(
  dmr_bed,
  file = "significant_DMRs_FDR05.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ------------------------------------------------------------
# 10. Print checks and results
# ------------------------------------------------------------

cat("\nINPUT CHECKS\n")
cat("Total WGBS regions:", nrow(dmr_all), "\n")
cat("Significant DMRs, FDR < 0.05:", length(dmr_gr), "\n")
cat("Kidney ChromHMM intervals:", length(kidney_chromhmm_gr), "\n")
cat("Liver ChromHMM intervals:", length(liver_chromhmm_gr), "\n")
cat("Kidney Het intervals:", length(kidney_het_gr), "\n")
cat("Liver Het intervals:", length(liver_het_gr), "\n")
cat("Kidney-specific ATAC peaks:", length(kidney_atac_gr), "\n")
cat("Liver-specific ATAC peaks:", length(liver_atac_gr), "\n")

cat("\nOVERLAP RESULTS\n")
print(summary_table, row.names = FALSE)

cat("\nOutput files created:\n")
cat("task4_overlap_summary.tsv\n")
cat("task4_DMR_overlap_details.csv\n")
cat("significant_DMRs_FDR05.bed\n")
