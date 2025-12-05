ECHOS Data Processing & Analysis Pipeline (v1.0)

This repository provides a modular and reproducible analysis pipeline for processing CUT&Tag / ECHOS datasets, performing peak calling, differential analysis, RNA–chromatin integration, and specialized micronuclei hotspot analyses.
Each module is designed to be stand-alone, allowing flexible reuse across histone marks, species, and experimental designs.

Pipeline Overview (Modules 01–07)
01_setup_and_basic_QC.*

Initial project setup and sample-level QC.

Functions:

Define project directory structure and sample metadata

Verify BAM integrity, optional read statistics

Output basic QC summaries

Inputs: aligned BAM files
Outputs: QC tables and optional plots

02_peak_calling_and_count_matrix.*

Genome-wide peak calling and peak × sample matrix construction.

Functions:

MACS2 peak calling (narrow / broad)

Construct a union peak set using bedtools merge

Generate peak counts (featureCounts or bedtools coverage)

Inputs: BAM files
Outputs:

*_peaks.narrowPeak or *_broadPeak

union peak BED

peak count matrix (CSV/TSV)

03_diff_peaks_edgeR_or_DESeq2.*

Differential peak analysis.

Functions:

Normalize peak counts

Perform differential testing (e.g., basal vs luminal)

Export log2FC, FDR, significance labels

Optional QC plots (MA plot, volcano, PCA)

Inputs: peak matrix
Outputs: differential peak tables + plots

04_homer_peak_annotation_and_motif.*

Peak annotation with HOMER.

Functions:

Annotate peaks relative to gene features

Link peaks to nearest genes

Optional: motif enrichment on selected peak subsets

Inputs: union peaks
Outputs: annotation tables (annotation.tsv) and optional motif folders

05a_diff_peaks_DESeq2.R

Peak-level DE using DESeq2 (recommended primary method).

Functions:

Construct DESeqDataSet

Fit model, export unshrunken and shrunken results

Generate MA/volcano/heatmap plots

Inputs: peak count matrix + sample metadata
Outputs:

deseq2_results_unshrunken.csv

optional shrunken results

QC plots (PDF)

06_rna_chromatin_integration.py

General RNA × chromatin integration module.

Functions:

Load RNA DE results (auto-detect log2FC/FDR columns)

Load ChIP/ECHOS differential peaks + annotation

Map peaks → genes (promoter-only or best peak)

Produce an integrated RNA–ChIP table

Compute quadrant categories (Up/Up, Down/Down, Opposite, NS)

Plot RNA vs ChIP scatter with marker labeling

Export quadrant gene lists

Inputs:

RNA DESeq2 outputs

ChIP peak DE + annotation

Outputs:

joined gene table (CSV)

quadrant gene lists (CSV)

RNA–ChIP scatter plot (PDF)

summary quadrant counts

07_micronuclei_chrY_hotspot_analysis.py

Dedicated pipeline for micronuclei-specific hotspot detection on chrY.

Functions:

chrY enrichment QC

Compute chrY coverage relative to autosomes

Barplot enrichment summary

Per-chromosome specificity

Compare normalized MN vs genomic signal

Export fold-change tables and diagnostic plots

SES-normalized log2(IP/CTL) processing

bamCompare with SES scaling and multi-scale smoothing

multiBigwigSummary matrices on chrY

Hotspot/domain detection

Compute R_MN, ΔR_TR, ΔR_UT

Threshold bins, merge contiguous domains, apply minimum length filters

Identify replicate-consistent MN hotspots

Visualization

Export track plots with highlighted MN-specific or replicate-consistent domains

Inputs: MN and whole-nucleus BAM files
Outputs:

chrY enrichment tables and plots

log2 matrices and ΔR tracks

hotspot BED files (5kb/50kb, 20kb/50kb versions)

SES log2 track PDF overlays