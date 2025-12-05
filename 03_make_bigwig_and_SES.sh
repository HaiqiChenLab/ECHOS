#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------------------
# 03_make_bigwig_and_SES.sh
#
# Universal signal-track generation script for CUT&Tag / ECHOS datasets.
#
# This script generates:
#   (1) Per-sample coverage bigWig from BAM (e.g. CPM-normalized)
#   (2) Signal tracks from MACS2 bedGraph outputs
#       - treat_pileup, control_lambda, and fold-enrichment (FE)
#   (3) Pairwise comparison tracks:
#       - log2(IP / control) from bigWig (bigwigCompare)
#       - SES-normalized log2(IP / control) from BAM (bamCompare + SES)
#   (4) Mean bigWig tracks by averaging multiple bigWigs (e.g., replicates)
#
# The script is marker- and project-agnostic; you only need to:
#   - Provide the BAM and MACS2 directories
#   - Provide a chrom.sizes file
#   - Edit the comparison / averaging arrays for your specific project
#
# Usage:
#   bash scripts/03_make_bigwig_and_SES.sh \
#       <bam_dir> \
#       <macs2_dir> \
#       <chrom_sizes> \
#       <out_base>
#
# Example:
#   bash scripts/03_make_bigwig_and_SES.sh \
#       results/align/bam \
#       results/macs2 \
#       config/hg38.chrom.sizes \
#       results/bigwig
#
# Output structure:
#   <out_base>/
#       from_bam/           # per-sample coverage bigWigs (e.g. CPM)
#       from_bam/ratio/     # pairwise log2 or SES-normalized tracks
#       from_peak/          # bigWigs derived from MACS2 bedGraphs
#       from_mean/          # mean tracks (averaged over multiple bigWigs)
#
# Requirements:
#   - deepTools >= 3.5 (bamCoverage, bigwigCompare, bamCompare, bigwigAverage)
#   - MACS2 (for bdgcmp to compute FE)
#   - bedGraphToBigWig (UCSC tools)
#   - bedClip (optional but recommended)
#   - samtools (only if you need to derive chrom.sizes from BAM; not required
#     if you already have a proper chrom.sizes file)
#
# --------------------------------------------------------------------

BAM_DIR="$1"
MACS2_DIR="$2"
CHROMSIZES="$3"
OUT_BASE="$4"

OUT_BAM_BW="${OUT_BASE}/from_bam"
OUT_BAM_RATIO="${OUT_BAM_BW}/ratio"
OUT_PEAK_BW="${OUT_BASE}/from_peak"
OUT_MEAN_BW="${OUT_BASE}/from_mean"

mkdir -p "$OUT_BAM_BW" "$OUT_BAM_RATIO" "$OUT_PEAK_BW" "$OUT_MEAN_BW"

# Optional: region and blacklist for SES (e.g., chrY-only micronucleus analysis)
# Set to empty strings ("") to disable.
SES_REGION="${SES_REGION:-}"        # e.g. "chrY"
SES_BLACKLIST="${SES_BLACKLIST:-}"  # e.g. "/path/to/blacklist.bed"

# --------------------------------------------------------------------
# Helper: load modules if available (for cluster environments)
# --------------------------------------------------------------------
load_project_modules() {
  if command -v module >/dev/null 2>&1; then
    module load deeptools/3.5.0 2>/dev/null || module load deepTools/3.5.0 2>/dev/null || true
    module load macs/2.1.2 2>/dev/null || true
    module load UCSC_userApps/v317 2>/dev/null || true
  fi
}
load_project_modules

# --------------------------------------------------------------------
# Helper: ensure command exists
# --------------------------------------------------------------------
ensure_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: command '$cmd' not found in PATH. Please install or load it." >&2
    exit 1
  fi
}

# Required tools
ensure_cmd bamCoverage
ensure_cmd bedGraphToBigWig
ensure_cmd macs2

# Optional tools (sections will be skipped if missing)
if ! command -v bigwigCompare >/dev/null 2>&1; then
  echo "[WARN] bigwigCompare not found. Pairwise bigWig log2 tracks will be skipped."
fi
if ! command -v bamCompare >/dev/null 2>&1; then
  echo "[WARN] bamCompare not found. SES-based tracks from BAM will be skipped."
fi
if ! command -v bigwigAverage >/dev/null 2>&1; then
  echo "[WARN] bigwigAverage not found. Mean bigWig generation will be skipped."
fi
if ! command -v bedClip >/dev/null 2>&1; then
  echo "[WARN] bedClip not found. bedGraph will not be clipped to chrom.sizes."
fi

# --------------------------------------------------------------------
# Helper: clip + sort bedGraph before conversion to bigWig
# --------------------------------------------------------------------
clip_and_sort_bdg() {
  local in_bdg="$1"
  local out_bdg_sorted="$2"
  local tmp_bdg="${out_bdg_sorted%.sorted.bdg}.tmpclip.bdg"

  if command -v bedClip >/dev/null 2>&1; then
    bedClip "$in_bdg" "$CHROMSIZES" "$tmp_bdg"
    LC_ALL=C sort -k1,1 -k2,2n "$tmp_bdg" > "$out_bdg_sorted"
    rm -f "$tmp_bdg"
  else
    # No clipping; sort only. bedGraphToBigWig may complain if coordinates exceed chrom.sizes.
    LC_ALL=C sort -k1,1 -k2,2n "$in_bdg" > "$out_bdg_sorted"
  fi
}

# --------------------------------------------------------------------
# SECTION 1: per-sample coverage bigWig from BAM (e.g., CPM)
# --------------------------------------------------------------------
echo "============================================================"
echo "[STEP 1] Generating coverage bigWig from BAM files (bamCoverage)"
echo "============================================================"

shopt -s nullglob
for bam in "${BAM_DIR}"/*.bam; do
  base=$(basename "$bam" .bam)
  out_bw="${OUT_BAM_BW}/${base}.cpm.bw"

  if [ -s "$out_bw" ]; then
    echo "[SKIP] Coverage bigWig already exists: $out_bw"
    continue
  fi

  echo "[INFO] bamCoverage -> $(basename "$out_bw")"
  bamCoverage \
    -b "$bam" \
    -o "$out_bw" \
    --normalizeUsing CPM \
    --binSize 25 \
    --centerReads \
    --extendReads \
    --ignoreDuplicates \
    --minMappingQuality 30
done
shopt -u nullglob

# --------------------------------------------------------------------
# SECTION 2: bigWig from MACS2 bedGraph (treat_pileup / control_lambda / FE)
# --------------------------------------------------------------------
echo
echo "============================================================"
echo "[STEP 2] Generating bigWig from MACS2 bedGraph outputs"
echo "============================================================"

if [ ! -s "$CHROMSIZES" ]; then
  echo "ERROR: chrom.sizes file not found or empty: $CHROMSIZES" >&2
  exit 1
fi

shopt -s nullglob
treat_files=( "${MACS2_DIR}"/*_treat_pileup.bdg )
if [ "${#treat_files[@]}" -eq 0 ]; then
  echo "[WARN] No *_treat_pileup.bdg files found in ${MACS2_DIR}. Skipping Section 2."
else
  for treat_bdg in "${treat_files[@]}"; do
    base=$(basename "$treat_bdg" _treat_pileup.bdg)
    ctrl_bdg="${MACS2_DIR}/${base}_control_lambda.bdg"

    t_sorted="${OUT_PEAK_BW}/${base}_treat_pileup.sorted.bdg"
    c_sorted="${OUT_PEAK_BW}/${base}_control_lambda.sorted.bdg"
    t_bw="${OUT_PEAK_BW}/${base}_treat_pileup.bw"
    c_bw="${OUT_PEAK_BW}/${base}_control_lambda.bw"

    echo "[INFO] Converting MACS2 pileup to bigWig: ${base}"

    # treat_pileup -> bigWig
    clip_and_sort_bdg "$treat_bdg" "$t_sorted"
    bedGraphToBigWig "$t_sorted" "$CHROMSIZES" "$t_bw"

    # control_lambda -> bigWig (if available)
    if [ -s "$ctrl_bdg" ]; then
      clip_and_sort_bdg "$ctrl_bdg" "$c_sorted"
      bedGraphToBigWig "$c_sorted" "$CHROMSIZES" "$c_bw"
    else
      echo "[WARN] control_lambda bedGraph not found for ${base}. FE computation will be skipped."
      continue
    fi

    # FE (fold-enrichment) -> bigWig
    fe_bdg="${OUT_PEAK_BW}/${base}_FE.bdg"
    fe_sorted="${OUT_PEAK_BW}/${base}_FE.sorted.bdg"
    fe_bw="${OUT_PEAK_BW}/${base}_FE.bw"

    echo "[INFO] Computing fold-enrichment (FE) and converting to bigWig: ${base}"
    macs2 bdgcmp -t "$treat_bdg" -c "$ctrl_bdg" -m FE -o "$fe_bdg"
    clip_and_sort_bdg "$fe_bdg" "$fe_sorted"
    bedGraphToBigWig "$fe_sorted" "$CHROMSIZES" "$fe_bw"
  done
fi
shopt -u nullglob

# --------------------------------------------------------------------
# SECTION 3A: log2(IP / control) from bigWig (bigwigCompare)
#
# This is what you used in the cervix H3K4me3 project:
#   - take CPM bigWigs from BAM
#   - compute log2(treatment / control) per pair
#
# You can define comparisons by editing the COMPARISONS_BW array below.
# Format of each entry:
#   "name|bw1|bw2"
#
# where:
#   name = output prefix (e.g. hCVX_K4_basal_1)
#   bw1  = bigWig file 1 (e.g. treatment CPM)
#   bw2  = bigWig file 2 (e.g. control CPM)
#
# Output:
#   <out_base>/from_bam/ratio/<name>.log2ratio.bw
# --------------------------------------------------------------------

echo
echo "============================================================"
echo "[STEP 3A] Pairwise log2(IP / control) from bigWig (bigwigCompare)"
echo "============================================================"

declare -a COMPARISONS_BW=(
  # Example (human cervix H3K4me3):
  # "hCVX_K4_basal_1|${OUT_BAM_BW}/hCVX_K4_basal_1_hg38.cpm.bw|${OUT_BAM_BW}/hCVX_who_basal_1_hg38.cpm.bw"
  # "hCVX_K4_basal_2|${OUT_BAM_BW}/hCVX_K4_basal_2_hg38.cpm.bw|${OUT_BAM_BW}/hCVX_who_basal_2_hg38.cpm.bw"
)

if ! command -v bigwigCompare >/dev/null 2>&1; then
  echo "[WARN] bigwigCompare not available. Skipping Section 3A."
else
  if [ "${#COMPARISONS_BW[@]}" -eq 0 ]; then
    echo "[INFO] No COMPARISONS_BW defined. Skipping Section 3A."
  else
    for item in "${COMPARISONS_BW[@]}"; do
      IFS='|' read -r name bw1 bw2 <<< "$item"

      if [ ! -s "$bw1" ] || [ ! -s "$bw2" ]; then
        echo "[WARN] Missing bigWig(s) for comparison '$name':"
        echo "       bw1 = $bw1"
        echo "       bw2 = $bw2"
        continue
      fi

      out_bw="${OUT_BAM_RATIO}/${name}.log2ratio.bw"
      echo "[INFO] bigwigCompare (log2) -> $(basename "$out_bw")"

      bigwigCompare \
        -b1 "$bw1" -b2 "$bw2" \
        --operation log2 \
        --pseudocount 1 \
        -o "$out_bw"
    done
  fi
fi

# --------------------------------------------------------------------
# SECTION 3B: SES-normalized log2(IP / control) from BAM (bamCompare)
#
# This captures your micronucleus analysis:
#   - SES-normalized, smoothed log2(IP / control), often restricted
#     to specific regions (e.g., chrY) and with optional blacklist.
#
# Define comparisons by editing the COMPARISONS_SES array.
# Format of each entry:
#   "name|bam1|bam2|binSize|smoothLength"
#
# where:
#   name         = output prefix (e.g. MN_K27ac_pss_1_vs_whole_pss_1)
#   bam1         = IP BAM
#   bam2         = control BAM
#   binSize      = bin size (e.g., 2000 or 5000)
#   smoothLength = smoothing length (e.g., 20000 or 25000)
#
# Optional globals:
#   SES_REGION    : if non-empty, added as `--region SES_REGION`
#   SES_BLACKLIST : if non-empty, added as `--blackListFileName SES_BLACKLIST`
#
# Output:
#   <out_base>/from_bam/ratio/<name>.SES.bw
# --------------------------------------------------------------------

echo
echo "============================================================"
echo "[STEP 3B] SES-normalized log2(IP / control) from BAM (bamCompare)"
echo "============================================================"

declare -a COMPARISONS_SES=(
  # Example (micronuclei H3K27ac on chrY, SES scaling):
  # "MN_K27ac_pss_1_vs_whole_pss_1|${BAM_DIR}/MN_K27ac_pss_1.bam|${BAM_DIR}/MN_whole_pss_1.bam|2000|20000"
  # "MN_K27ac_pss_2_vs_whole_pss_2|${BAM_DIR}/MN_K27ac_pss_2.bam|${BAM_DIR}/MN_whole_pss_2.bam|2000|20000"
)

if ! command -v bamCompare >/dev/null 2>&1; then
  echo "[WARN] bamCompare not available. Skipping Section 3B."
else
  if [ "${#COMPARISONS_SES[@]}" -eq 0 ]; then
    echo "[INFO] No COMPARISONS_SES defined. Skipping Section 3B."
  else
    # Build optional region/blacklist args
    ses_region_arg=""
    ses_black_arg=""
    if [ -n "$SES_REGION" ]; then
      ses_region_arg="--region ${SES_REGION}"
      echo "[INFO] SES region: ${SES_REGION}"
    fi
    if [ -n "$SES_BLACKLIST" ]; then
      ses_black_arg="--blackListFileName ${SES_BLACKLIST}"
      echo "[INFO] SES blacklist: ${SES_BLACKLIST}"
    fi

    for item in "${COMPARISONS_SES[@]}"; do
      IFS='|' read -r name bam1 bam2 binSize smoothLen <<< "$item"

      if [ ! -s "$bam1" ] || [ ! -s "$bam2" ]; then
        echo "[WARN] Missing BAM(s) for SES comparison '$name':"
        echo "       bam1 = $bam1"
        echo "       bam2 = $bam2"
        continue
      fi

      out_bw="${OUT_BAM_RATIO}/${name}.SES.bw"
      echo "[INFO] bamCompare (log2 + SES) -> $(basename "$out_bw")"

      bamCompare \
        -b1 "$bam1" -b2 "$bam2" \
        --operation log2 \
        --pseudocount 1 \
        --scaleFactorsMethod SES \
        --binSize "$binSize" \
        --smoothLength "$smoothLen" \
        --minMappingQuality 30 \
        --extendReads 150 \
        $ses_region_arg \
        $ses_black_arg \
        -o "$out_bw"
    done
  fi
fi

# --------------------------------------------------------------------
# SECTION 4: mean bigWig tracks (averaged replicates)
#
# Uses deepTools "bigwigAverage" to compute the mean track across
# multiple bigWig files. This is analogous to your barr body averaging
# script that used bigwigMerge + averaging, but here we use the
# higher-level deepTools interface.
#
# Define groups by editing AVERAGE_GROUPS:
#   "out_name|bw1,bw2[,bw3,...]"
#
# Example (barr body H3K27me3, RPGC bigWigs):
#   "Young_barr_RPGC50bp_mean|${BW_DIR}/Young_barr_1_S1_hg38_RPGC_50bp.bw,${BW_DIR}/Young_barr_2_S2_hg38_RPGC_50bp.bw"
# --------------------------------------------------------------------

echo
echo "============================================================"
echo "[STEP 4] Mean bigWig tracks (averaged replicates)"
echo "============================================================"

declare -a AVERAGE_GROUPS=(
  # Example (barr body H3K27me3):
  # "Young_barr_RPGC50bp_mean|${OUT_BAM_BW}/Young_barr_1_S1_hg38_RPGC_50bp.bw,${OUT_BAM_BW}/Young_barr_2_S2_hg38_RPGC_50bp.bw"
  # "Old_barr_RPGC50bp_mean|${OUT_BAM_BW}/Old_barr_1_S6_hg38_RPGC_50bp.bw,${OUT_BAM_BW}/Old_barr_2_S7_hg38_RPGC_50bp.bw"
)

if ! command -v bigwigAverage >/dev/null 2>&1; then
  echo "[WARN] bigwigAverage not available. Skipping mean bigWig generation."
else
  if [ "${#AVERAGE_GROUPS[@]}" -eq 0 ]; then
    echo "[INFO] No AVERAGE_GROUPS defined. Skipping Section 4."
  else
    for item in "${AVERAGE_GROUPS[@]}"; do
      IFS='|' read -r out_name bw_list <<< "$item"
      IFS=',' read -r -a bws <<< "$bw_list"

      out_bw="${OUT_MEAN_BW}/${out_name}.bw"
      echo "[INFO] Generating mean bigWig: $(basename "$out_bw")"
      echo "       Inputs:"

      for bw in "${bws[@]}"; do
        echo "         - $bw"
        if [ ! -s "$bw" ]; then
          echo "ERROR: missing or empty bigWig in mean group '$out_name': $bw" >&2
          exit 1
        fi
      done

      bigwigAverage \
        -b "${bws[@]}" \
        -o "$out_bw"
    done
  fi
fi

echo
echo "[DONE] All bigWig/SES generation steps completed."
echo "Output base directory: $OUT_BASE"
