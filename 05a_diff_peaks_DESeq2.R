#!/usr/bin/env Rscript

########################################################################
# 05a_diff_peaks_DESeq2.R
#
# Generic differential peak analysis with DESeq2.
#
# This script:
#   1) Reads a peak-by-sample count matrix and sample metadata.
#   2) Builds a DESeq2 model with a simple design (~ condition).
#   3) Performs a specified 2-group contrast (contrast_level vs ref_level).
#   4) Applies log2 fold-change shrinkage (apeglm if available).
#   5) Writes result tables and basic QC plots (MA, PCA, volcano).
#
# Intended to be used after:
#   05_diff_peaks_prepare_counts.py
#
# Edit ONLY the CONFIG section below for each contrast.
########################################################################

## =========================
## CONFIG (EDIT THIS BLOCK)
## =========================

# Path to counts matrix (rows = peaks, columns = samples)
# Produced by 05_diff_peaks_prepare_counts.py
counts_file <- "results/diff_peaks/example_contrast/counts_peak_matrix.csv"

# Path to sample metadata table (one row per sample)
# Must contain at least: sample, group (or another condition column)
meta_file   <- "results/diff_peaks/example_contrast/sample_metadata.csv"

# Output base directory and project tag (for subfolder naming)
outdir_base <- "results/diff_peaks"
project_tag <- "example_contrast"   # e.g. "H3K4me3_cervix_basal_vs_luminal"

# Column in metadata used as condition (e.g. "group", "condition")
condition_col <- "group"

# Reference and contrast levels (factor levels of condition_col)
# Example 1: basal vs luminal
#   ref_level      <- "basal"
#   contrast_level <- "luminal"
#
# Example 2: young vs aged
#   ref_level      <- "young"
#   contrast_level <- "aged"
ref_level      <- "basal"
contrast_level <- "luminal"

# Filtering thresholds:
# Keep peaks with >= min_counts in at least min_samples samples.
min_counts  <- 10
min_samples <- 2

# Prefix used in output filenames
result_prefix <- "deseq2_peaks"

## =========================
## END OF CONFIG
## =========================


suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(grDevices)
})

has_apeglm <- suppressWarnings(requireNamespace("apeglm", quietly = TRUE))

# Output directories
outdir  <- file.path(outdir_base, project_tag)
resdir  <- file.path(outdir, "deseq2_results")
plotdir <- file.path(outdir, "plots")
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "]", ..., "\n")

say("Project tag:", project_tag)
say("Counts file:", counts_file)
say("Meta file  :", meta_file)

# ----------------------
# 1) Read input tables
# ----------------------

if (!file.exists(counts_file)) stop("Counts matrix not found: ", counts_file)
if (!file.exists(meta_file))   stop("Sample metadata not found: ", meta_file)

counts_raw <- read.csv(counts_file, row.names = 1, check.names = FALSE)
meta       <- read.csv(meta_file, stringsAsFactors = FALSE)

if (!all(c("sample", condition_col) %in% colnames(meta))) {
  stop(sprintf(
    "Metadata must contain columns: 'sample' and '%s'. Found: %s",
    condition_col, paste(colnames(meta), collapse = ", ")
  ))
}

# Clean column names of counts (light normalization; safe for HOMER/deepTools outputs)
cn <- colnames(counts_raw)
cn <- sub("^.*/", "", cn)                                   # drop path prefixes
cn <- sub("[[:space:]]*Tag Count in given bp.*$", "", cn)   # drop HOMER suffix
cn <- sub("_hg[0-9]+$", "", cn)                             # drop trailing _hg19/_hg38 if present
cn <- trimws(cn)
colnames(counts_raw) <- cn

# Clean sample names in metadata
meta$sample <- trimws(sub("_hg[0-9]+$", "", meta$sample))

# Keep overlapping samples and align order
common <- intersect(colnames(counts_raw), meta$sample)
if (length(common) < 2) {
  message("[DEBUG] counts columns:")
  print(colnames(counts_raw))
  message("[DEBUG] meta$sample:")
  print(meta$sample)
  stop("< 2 overlapping samples between counts and metadata")
}

counts <- counts_raw[, common, drop = FALSE]
meta   <- meta[match(colnames(counts), meta$sample), , drop = FALSE]

say("Number of samples:", ncol(counts))
say("Samples:", paste(colnames(counts), collapse = ", "))

# ------------------------------
# 2) Prepare condition factor
# ------------------------------

if (!condition_col %in% colnames(meta)) {
  stop("Condition column not found in metadata: ", condition_col)
}

cond <- factor(meta[[condition_col]])
if (!all(c(ref_level, contrast_level) %in% levels(cond))) {
  stop(sprintf(
    "Requested levels not found in '%s': ref_level='%s', contrast_level='%s'. Levels: %s",
    condition_col, ref_level, contrast_level, paste(levels(cond), collapse = ", ")
  ))
}

cond <- relevel(cond, ref = ref_level)
meta[[condition_col]] <- cond

say("Condition column:", condition_col)
say("Levels:", paste(levels(cond), collapse = ", "))
say("Reference level:", ref_level, "Contrast level:", contrast_level)

# ------------------------------
# 3) Convert counts to integers
# ------------------------------

counts_mat <- as.matrix(counts)
storage.mode(counts_mat) <- "numeric"
counts_mat <- round(counts_mat)
storage.mode(counts_mat) <- "integer"

# ------------------------------
# 4) Basic filtering (low counts)
# ------------------------------

keep <- rowSums(counts_mat >= min_counts) >= min_samples
counts_f <- counts_mat[keep, , drop = FALSE]
say(sprintf("Kept peaks after filtering: %d / %d", nrow(counts_f), nrow(counts_mat)))

if (nrow(counts_f) < 5) {
  warning("Fewer than 5 peaks left after filtering; results may not be meaningful.")
}

# ------------------------------
# 5) Build DESeq2 object & run
# ------------------------------

coldata <- data.frame(
  row.names  = meta$sample,
  condition  = meta[[condition_col]]
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_f,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)

# Unshrunken results
res_raw <- results(dds, contrast = c("condition", contrast_level, ref_level))

# LFC shrinkage
coefName <- grep("condition_", resultsNames(dds), value = TRUE)[1]
if (is.na(coefName)) {
  # Fallback: first term containing 'condition'
  coefName <- resultsNames(dds)[grep("condition", resultsNames(dds))[1]]
}

if (has_apeglm && !is.na(coefName)) {
  say("Using apeglm for LFC shrinkage (coef =", coefName, ") ...")
  res_shrunk <- tryCatch(
    lfcShrink(dds, coef = coefName, type = "apeglm"),
    error = function(e) {
      say("apeglm shrinkage failed:", conditionMessage(e), "-> fallback to 'normal'")
      tryCatch(
        lfcShrink(dds, coef = coefName, type = "normal"),
        error = function(e2) res_raw
      )
    }
  )
} else {
  say("apeglm not available; using 'normal' shrinkage (or raw LFC if that fails)...")
  res_shrunk <- tryCatch(
    lfcShrink(dds, coef = coefName, type = "normal"),
    error = function(e) res_raw
  )
}

res_df <- as.data.frame(res_shrunk)
res_df$peak_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

res_raw_df <- as.data.frame(res_raw)
res_raw_df$peak_id <- rownames(res_raw_df)
res_raw_df <- res_raw_df[order(res_raw_df$padj), ]

# ------------------------------
# 6) Write result tables
# ------------------------------

say("Writing result tables ...")

res_shrunk_file <- file.path(
  resdir,
  sprintf("%s_%s_vs_%s_shrunken.csv", result_prefix, contrast_level, ref_level)
)
res_raw_file <- file.path(
  resdir,
  sprintf("%s_%s_vs_%s_unshrunken.csv", result_prefix, contrast_level, ref_level)
)
norm_counts_file <- file.path(
  resdir,
  sprintf("%s_normalized_counts.csv", result_prefix)
)

write.csv(res_df, res_shrunk_file, row.names = FALSE)

raw_out <- res_raw_df
colnames(raw_out) <- sub("^log2FoldChange$", "log2FoldChange_raw", colnames(raw_out))
colnames(raw_out) <- sub("^padj$", "padj_raw", colnames(raw_out))
write.csv(raw_out, res_raw_file, row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, norm_counts_file)

# rlog / vst (optional, may fail for very large objects)
rlog_file <- file.path(resdir, sprintf("%s_rlog_matrix.csv", result_prefix))
vst_file  <- file.path(resdir, sprintf("%s_vst_matrix.csv", result_prefix))

try({
  rld <- rlog(dds, blind = FALSE)
  write.csv(assay(rld), rlog_file)
}, silent = TRUE)

try({
  vstd <- vst(dds, blind = FALSE)
  write.csv(assay(vstd), vst_file)
}, silent = TRUE)

# ------------------------------
# 7) Plots (MA, PCA, Volcano)
# ------------------------------

# --- MA plot (using shrunken LFC) ---
ma_pdf <- file.path(
  plotdir,
  sprintf("MA_%s_vs_%s.pdf", contrast_level, ref_level)
)
pdf(ma_pdf, width = 6, height = 4)
plotMA(
  res_shrunk,
  ylim = c(-4, 4),
  main = sprintf("MA plot (%s vs %s)", contrast_level, ref_level)
)
dev.off()

# --- PCA plot (rlog) ---
pca_pdf <- file.path(
  plotdir,
  sprintf("PCA_rlog_%s_vs_%s.pdf", contrast_level, ref_level)
)
try({
  if (!exists("rld")) rld <- rlog(dds, blind = TRUE)
  p <- plotPCA(rld, intgroup = "condition") +
    ggtitle(sprintf("PCA (rlog) — %s vs %s", contrast_level, ref_level))
  ggplot2::ggsave(filename = pca_pdf, plot = p,
                  width = 6, height = 4, units = "in")
}, silent = TRUE)

# --- Volcano plot (shrunken LFC + padj) ---
vol <- na.omit(as.data.frame(res_shrunk))
vol$signif <- with(
  vol,
  ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1,
         "Significant", "NS")
)

volcano_pdf <- file.path(
  plotdir,
  sprintf("Volcano_%s_vs_%s.pdf", contrast_level, ref_level)
)

pvol <- ggplot(vol, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = c("Significant" = "#D62728", "NS" = "#7F7F7F")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(
    title = sprintf("Volcano: %s vs %s", contrast_level, ref_level),
    x     = sprintf("log2FC (%s / %s)", contrast_level, ref_level),
    y     = "-log10(FDR)",
    color = NULL
  ) +
  theme_minimal(base_size = 11)

ggplot2::ggsave(volcano_pdf, plot = pvol, width = 6, height = 4, units = "in")

# ------------------------------
# 8) Session info
# ------------------------------

session_file <- file.path(
  resdir,
  sprintf("%s_sessionInfo.txt", result_prefix)
)
writeLines(capture.output(sessionInfo()), con = session_file)

say("DONE.")
say("Results dir:", resdir)
say("Plots   dir:", plotdir)
