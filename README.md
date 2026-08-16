# ATAC-seq Analysis of Fetal Mouse Kidney and Liver

End-to-end chromatin accessibility analysis covering the nf-core/atacseq pipeline through integration with DNA methylation, histone modifications and gene expression, applied to eight ATAC-seq libraries from developing mouse tissue (Gorkin et al., 2020).

## Overview

This project processes eight paired-end ATAC-seq libraries from mouse fetal **liver** and **kidney** at embryonic days **E14.5** and **E15.5** (two biological replicates per tissue and timepoint) to identify tissue-specific regulatory elements, infer the transcription factors driving them, and relate accessibility to the other epigenetic layers profiled in parallel by the other course groups. Data are restricted to a reduced mm10 reference containing chromosomes 18 and 19 only (145,658,490 bp).

Completed as part of Computational Methods for Epigenome Analysis (COMPEPIWS 2026) at Saarland University, Walter Lab & Müller Lab (ICBB), group `atacseq1`.

**Authors:** Divyashree Mohan and Ajay Nimbalkar · Mentor: Midhuna Immaculate Maran

**Samples:** `kidney_14_5_REP1/2`, `kidney_15_5_REP1/2`, `liver_14_5_REP1/2`, `liver_15_5_REP1/2` — 1.7–4.4 M reads each

## Pipeline

1. **Read processing & alignment** — nf-core/atacseq v2.1.2 with Singularity: cutadapt trimming, BWA-MEM alignment, Picard MarkDuplicates, ENCODE mm10 blacklist filtering
2. **Peak calling** — MACS2 narrow-peak mode, 1,456–3,415 peaks per replicate merged into a 6,259-peak consensus set
3. **Quality control** — FastQC and MultiQC; FRiP, TSS enrichment, fragment-size periodicity and promoter-overlap metrics
4. **Signal summarisation** — ChrAccR over three region sets (genome-wide tiling, promoters, consensus peaks)
5. **Filtering & normalisation** — low-coverage and chrM removal, then quantile normalisation
6. **Exploratory analysis** — PCA and UMAP per region set, clustered heatmaps of the most variable regions
7. **Motif activity** — chromVAR with 579 JASPAR vertebrate motifs, variability ranking and deviation Z-scores
8. **Differential accessibility** — kidney vs liver on raw counts (padj < 0.01, |log2FC| > 3), volcano plot, nearest-TSS annotation, Wilcoxon test on TSS distances
9. **TF footprinting** — aggregate Tn5 insertion profiles at GATA1, GATA4, GBX2 and HOXD13 motifs from tissue-merged datasets
10. **Locus-level integration** — IGV session combining accessibility, ChromHMM segmentation, WGBS signal and RNA-seq
11. **Quantitative integration** — methylation per chromatin state; methylation vs accessibility across 5,720 exactly matched peaks; DMR overlap with heterochromatin and differential peaks; DMR-to-gene assignment vs expression

## Key results

- **Tissue dominates over developmental time** — kidney and liver separate cleanly on consensus peaks (PC1 = 79.69% of variance), while E14.5 and E15.5 samples intermix freely within each tissue (`figures/chraccr/08_pca_consensus_peaks.png`)
- **Differential accessibility is distal** — 509 peaks more accessible in liver vs 109 in kidney, but almost none at promoters (7 vs 2); kidney-specific peaks sit significantly farther from the nearest TSS than non-differential peaks (Wilcoxon p = 2.32 × 10⁻¹⁰) (`figures/chraccr/19_volcano_differential_accessibility.png`, `20_violin_tss_distance.png`)
- **One TF family dominates motif variability** — the five most variable motifs are all GATA-family (Gata4 9.18, Gata1 7.80, GATA6 7.26, GATA3 6.72, GATA1::TAL1 6.47), reflecting both hepatocyte identity (GATA4/6) and the erythroid progenitor population of the fetal liver (`figures/chraccr/17_motif_variability.png`)
- **Motif enrichment is not occupancy** — footprinting confirms tissue-specific binding at GATA sites, but the kidney-enriched HOXD13 and GBX2 motifs show identical aggregate profiles in both tissues (`figures/chraccr/21_footprint_gata1.png`, `22_footprint_hoxd13.png`)
- **Methylation and accessibility are inversely related** — r = −0.748 (kidney) and −0.523 (liver) across 5,720 matched peaks; restricting to regions differential in both assays tightens this to 535 of 538 overlaps (99.44%) showing the inverse pattern (`figures/integrative/38_scatter_kidney_meth_vs_acc.png`, `41_dmr_atac_overlap.png`)
- **The methylation–expression relationship is context-dependent** — active promoters are unmethylated in both tissues (Pr_A: 0.009 / 0.014) while transcribed gene bodies are the most methylated compartment (Tx_S: 0.898 / 0.812), so gene-averaged methylation correlates only weakly with expression (ρ = −0.10 to −0.12) (`figures/integrative/36_methylation_by_chromatin_state.png`)

## Repository structure

```
report/              analysis_report.md    — full write-up: methods, 13 figures, interpretation
protocol/            part1.md, part2.md    — complete task-by-task workshop notebook
ChrAccR_analysis/    Browsable ChrAccR HTML report (summary, filtering,
                     normalization, exploratory, differential)
nextflow_run/        samplesheet.csv, config, MACS2 narrowPeak/summits/consensus
integrative/         task2–task5 R scripts, logs and result tables
figures/             36 figures in pipeline_qc/, chraccr/, integrative/
```

**[View the ChrAccR report](https://divyashreemohan05.github.io/ATACseq-Epigenomics-SS26/ChrAccR_analysis/index.html)** · **[Read the analysis report](report/analysis_report.md)**

## Notes on the data

Alignment files (BAM, bigWig) and MACS2 pileup bedGraphs are excluded — they run to several GB and are regenerable from the pipeline. The 200 bp aggregate count tracks used for IGV can be recreated with `exportCountTracks()` as documented in `protocol/part1.md`. Two ChrAccR outputs are also omitted for size: `diffObj_tiling.rds` (179 MB) and `diffTab_1_tiling.tsv` (59 MB).

WGBS, ChIP-seq and RNA-seq inputs used in the integrative sections were produced by the `wgbs1`, `wgbs2`, `chipseq1` and `rnaseq1` groups and shared through the course volume.

## Tools

nf-core/atacseq 2.1.2, Nextflow, Singularity, BWA-MEM, Picard, MACS2, MultiQC, FastQC · R with ChrAccR, chromVAR, DESeq2, GenomicRanges, edgeR, ComplexHeatmap, pheatmap, ggplot2, dplyr · bedtools, deepTools, IGV, ChromHMM · SLURM on the de.NBI cloud

## References

Gorkin, D. U. et al. (2020). The dynamic landscape of open chromatin during mouse embryonic development. *Nature* 583, 744–751.

Schep, A. N. et al. (2017). chromVAR: inferring transcription-factor-associated accessibility from single-cell epigenomic data. *Nature Methods* 14, 975–978.

[nf-core/atacseq](https://nf-co.re/atacseq) · [ChrAccR](https://epigenomeinformatics.github.io/ChrAccR/articles/overview.html)
