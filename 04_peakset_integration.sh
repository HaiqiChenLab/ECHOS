#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------------------
# 04_peakset_integration.sh
#
# Peak-set integration & basic QC metrics for CUT&Tag / ECHOS projects.
#
# This script is marker- and project-agnostic. It provides:
#   (1) Peak merging and union/intersection across methods / replicates
#       - bedtools sort + merge
#       - bedtools multiinter (multiIntersectBed)
#   (2) FRiP calculation for PE and SE datasets
#   (3) Signal matrix via deepTools multiBigwigSummary
#       (for correlation / heatmap / PCA plots)
#   (4) Placeholder hooks for NSC/RSC (phantompeakqualtools or ENCODE QC)
#
# The script is intentionally modular: you only need to edit the
# configuration section for your project (paths and sample lists).
#
# Usage:
#   bash scripts/04_peakset_integration.sh \
#       <out_dir>
#
# Example:
#   bash scripts/04_peakset_integration.sh results/peak_qc
#
# Downstream plotting (scatter, Venn, heatmaps, etc.) is handled by
# separate R/Python scripts (e.g., 04a_QC_plots_T7vsPSS.py).
#
# --------------------------------------------------------------------

OUT_DIR="$1"
mkdir -p "$OUT_DIR"

# --------------------------------------------------------------------
# 0. Load modules / ensure tools (adapt for your cluster)
# --------------------------------------------------------------------
if command -v module >/dev/null 2>&1; then
  module load bedtools2/2.31.1 2>/dev/null || module load bedtools/2.29.2 2>/dev/null || true
  module load samtools/1.19.2 2>/dev/null || module load samtools 2>/dev/null || true
  module load deeptools/3.5.0 2>/dev/null || module load deepTools/3.5.0 2>/dev/null || true
fi

ensure_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: command '$cmd' not found. Please install or module load it." >&2
    exit 1
  fi
}

ensure_cmd bedtools
ensure_cmd samtools
ensure_cmd multiBigwigSummary

# --------------------------------------------------------------------
# 1. Project-specific configuration (EDIT THIS SECTION)
# --------------------------------------------------------------------

# 1.1 Peak BED files for each method (example: H3K4me3 in HeLa)
# You can adapt these for other histone marks or species.

# Example: merged/intersected peaks per method (single BED per method)
PEAKS_METHOD_A=""  # e.g., ECHOS / PSS merged peaks
PEAKS_METHOD_B=""  # e.g., ChIP merged peaks
PEAKS_METHOD_C=""  # e.g., CUT&Tag merged peaks

# Example: replicate-level peak BEDs per method (for reproducibility)
# Format: "label:/path/to/rep1.bed,/path/to/rep2.bed"
REPLICATE_PEAK_SETS=(
  # "ECHOS:/path/to/echos_rep1.bed,/path/to/echos_rep2.bed"
  # "ChIP:/path/to/chip_rep1.bed,/path/to/chip_rep2.bed"
  # "CUTTag:/path/to/cutntag_rep1.bed,/path/to/cutntag_rep2.bed"
)

# 1.2 FRiP inputs: BAM + corresponding peak set, and whether the tech is PE or SE
# Format: "sample_id|bam|peaks|mode"
#   mode = PE or SE
FRIP_SAMPLES=(
  # Example (adapting from your FRiP script):
  # "CUTTag_rep1|/path/to/cutntag_rep1.sorted.bam|/path/to/CUTTag_peaks.bed|PE"
  # "PSS_rep1|/path/to/sorted_PSS_rep1.bam|/path/to/PSS_peaks.bed|PE"
  # "ChIP_rep1|/path/to/ChIP_rep1.sorted.bam|/path/to/ChIP_peaks.bed|SE"
)

# 1.3 multiBigwigSummary inputs: list of bigWigs and output prefix
BIGWIG_LIST=(
  # e.g.,
  # "/path/to/ChIP_rep1_treat_pileup.bw"
  # "/path/to/ChIP_rep2_treat_pileup.bw"
  # "/path/to/ECHOS_rep1_treat_pileup.bw"
  # "/path/to/ECHOS_rep2_treat_pileup.bw"
)
MULTIBW_BIN_SIZE=10000
MULTIBW_OUT_PREFIX="${OUT_DIR}/signal_matrix"  # will produce .npz and .tab

# --------------------------------------------------------------------
# 2. Peak merging & reproducibility (union / multiinter)
# --------------------------------------------------------------------
echo "============================================================"
echo "[STEP 2] Peak merging & reproducibility (bedtools)"
echo "============================================================"

# 2.1 Merge replicate peaks per method (if REPLICATE_PEAK_SETS defined)
MERGED_DIR="${OUT_DIR}/merged_peaks"
mkdir -p "$MERGED_DIR"

for item in "${REPLICATE_PEAK_SETS[@]}"; do
  IFS=':' read -r label bed_list <<< "$item"
  IFS=',' read -r -a beds <<< "$bed_list"

  out_merged="${MERGED_DIR}/${label}_replicates_merged.bed"
  out_union="${MERGED_DIR}/${label}_replicates_union.bed"

  echo "[INFO] Merging replicate peaks for: ${label}"
  cat "${beds[@]}" \
    | bedtools sort -i - \
    | bedtools merge -i - \
    > "$out_merged"

  # For convenience, treat merged = union here; if you need more complex logic,
  # you can modify this section.
  cp "$out_merged" "$out_union"
done

# 2.2 Multi-intersect across methods (if PEAKS_METHOD_* are set)
MULTIINTER_OUT="${OUT_DIR}/peak_multiinter.tsv"
if [[ -n "${PEAKS_METHOD_A}" || -n "${PEAKS_METHOD_B}" || -n "${PEAKS_METHOD_C}" ]]; then
  echo "[INFO] Running bedtools multiinter across methods..."
  # Build argument list dynamically
  args=()
  label_line="#chrom\tstart\tend"
  idx=1
  for peaks in "$PEAKS_METHOD_A" "$PEAKS_METHOD_B" "$PEAKS_METHOD_C"; do
    if [[ -n "$peaks" ]]; then
      args+=( -i "$peaks" )
      label_line="${label_line}\tmethod${idx}_count"
      ((idx++))
    fi
  done

  if (( ${#args[@]} > 0 )); then
    {
      echo -e "$label_line"
      bedtools multiinter "${args[@]}" \
        | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}'
    } > "$MULTIINTER_OUT"
    echo "[INFO] multiinter result: $MULTIINTER_OUT"
  fi
else
  echo "[INFO] No PEAKS_METHOD_* defined. Skipping multiinter step."
fi

# --------------------------------------------------------------------
# 3. FRiP calculation (PE vs SE-aware)
# --------------------------------------------------------------------
echo
echo "============================================================"
echo "[STEP 3] FRiP calculation"
echo "============================================================"

FRIP_RAW="${OUT_DIR}/FRiP_raw_counts.tsv"
FRIP_SUMMARY="${OUT_DIR}/FRiP_summary.tsv"
echo -e "sample\ttotal_frag_or_read\tin_peaks" > "$FRIP_RAW"

frip_pe() {
  local BAM="$1"
  local PEAKS="$2"
  local SAMPLE="$3"

  if [ ! -f "$BAM" ]; then
    echo "[SKIP] BAM not found: $BAM" >&2
    return 0
  fi
  if [ ! -f "$PEAKS" ]; then
    echo "[SKIP] Peaks not found: $PEAKS" >&2
    return 0
  fi

  # PE: count fragments via R1 (0x40)
  local TOTAL
  local IN_PEAKS
  TOTAL=$(samtools view -c -f 0x40 -F 4 "$BAM")
  IN_PEAKS=$(bedtools intersect -u -abam "$BAM" -b "$PEAKS" \
             | samtools view -c -f 0x40 -F 4)

  echo -e "${SAMPLE}\t${TOTAL}\t${IN_PEAKS}" >> "$FRIP_RAW"
}

frip_se() {
  local BAM="$1"
  local PEAKS="$2"
  local SAMPLE="$3"

  if [ ! -f "$BAM" ]; then
    echo "[SKIP] BAM not found: $BAM" >&2
    return 0
  fi
  if [ ! -f "$PEAKS" ]; then
    echo "[SKIP] Peaks not found: $PEAKS" >&2
    return 0
  fi

  # SE: count all mapped reads
  local TOTAL
  local IN_PEAKS
  TOTAL=$(samtools view -c -F 4 "$BAM")
  IN_PEAKS=$(bedtools intersect -u -abam "$BAM" -b "$PEAKS" \
             | samtools view -c -F 4)

  echo -e "${SAMPLE}\t${TOTAL}\t${IN_PEAKS}" >> "$FRIP_RAW"
}

# Ensure BAM indices
for item in "${FRIP_SAMPLES[@]}"; do
  IFS='|' read -r sid bam peaks mode <<< "$item"
  if [ -f "$bam" ] && [ ! -f "${bam}.bai" ]; then
    echo "[Index] $bam"
    samtools index "$bam"
  fi
done

# Run FRiP for all samples
for item in "${FRIP_SAMPLES[@]}"; do
  IFS='|' read -r sid bam peaks mode <<< "$item"
  echo "[FRiP] ${sid} (${mode})"
  case "$mode" in
    PE) frip_pe "$bam" "$peaks" "$sid" ;;
    SE) frip_se "$bam" "$peaks" "$sid" ;;
    *)  echo "[WARN] Unknown mode for FRiP: $mode (sample $sid)";;
  esac
done

# Add FRiP fraction
{
  echo -e "sample\ttotal_frag_or_read\tin_peaks\tFRiP"
  awk 'NR>1{frip=($2>0)?$3/$2:0; print $1"\t"$2"\t"$3"\t"frip}' "$FRIP_RAW"
} > "$FRIP_SUMMARY"
echo "[INFO] FRiP summary written to: $FRIP_SUMMARY"

# --------------------------------------------------------------------
# 4. multiBigwigSummary: signal matrix for correlation / heatmap / PCA
# --------------------------------------------------------------------
echo
echo "============================================================"
echo "[STEP 4] multiBigwigSummary (signal matrix)"
echo "============================================================"

if [ "${#BIGWIG_LIST[@]}" -eq 0 ]; then
  echo "[INFO] No BIGWIG_LIST defined. Skipping multiBigwigSummary."
else
  BW_ARGS=()
  for bw in "${BIGWIG_LIST[@]}"; do
    if [ -s "$bw" ]; then
      BW_ARGS+=( -b "$bw" )
    else
      echo "[WARN] bigWig not found or empty: $bw" >&2
    fi
  done

  if (( ${#BW_ARGS[@]} == 0 )); then
    echo "[WARN] No valid bigWigs found for multiBigwigSummary. Skipping."
  else
    NPZ="${MULTIBW_OUT_PREFIX}.npz"
    TAB="${MULTIBW_OUT_PREFIX}.tab"
    echo "[INFO] multiBigwigSummary -> ${NPZ}, ${TAB}"
    multiBigwigSummary bins \
      "${BW_ARGS[@]}" \
      --binSize "$MULTIBW_BIN_SIZE" \
      --outRawCounts "$TAB" \
      -out "$NPZ"
  fi
fi

# --------------------------------------------------------------------
# 5. Placeholder: NSC/RSC metrics (phantompeakqualtools / ENCODE)
# --------------------------------------------------------------------
echo
echo "============================================================"
echo "[STEP 5] NSC/RSC (placeholder)"
echo "============================================================"
echo "NOTE: NSC/RSC calculation is project-specific and usually depends"
echo "on phantompeakqualtools or the ENCODE ChIP-seq QC pipeline."
echo "You can add your own wrapper here to compute NSC/RSC for each BAM"
echo "and write the results to:"
echo "  ${OUT_DIR}/NSC_RSC_metrics.tsv"

echo
echo "[DONE] Peak-set integration & QC finished."
echo "Output directory: $OUT_DIR"
