#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------------------
# 02_macs2_callpeaks.sh
#
# Universal MACS2 peak-calling script for CUT&Tag / ECHOS datasets.
#
# This script:
#   • Reads sample metadata (sample → control relationships) from a TSV file
#   • Reads histone-mark–specific MACS2 parameters from a separate TSV file
#   • Calls MACS2 peaks for all “treatment” samples in a standardized manner
#
# The script is fully general:
#   • Works for any histone mark (H3K27ac, H3K4me3, H3K27me3, etc.)
#   • Works for any genome build (hg19, hg38, T2T-CHM13, mm10, ...)
#   • Works for any number of samples
#
# No assumptions are made about:
#   • experimental design
#   • sample naming
#   • control selection strategy
#
# Inputs:
#   1) samples.tsv
#   2) peakcalling_params.tsv
#   3) directory containing BAM files (from Step 01)
#   4) output directory for MACS2 results
#
# Usage example:
#   bash scripts/02_macs2_callpeaks.sh \
#       config/samples.tsv \
#       config/peakcalling_params.tsv \
#       results/align/bam \
#       results/macs2
#
# Expected samples.tsv (tab-delimited):
#   sample_id    fastq1   fastq2   mark     group       control_sample
#   H3K27ac_1    ...      ...      H3K27ac  treatment   Tn5_ctrl_1
#   Tn5_ctrl_1   ...      ...      control  control     NA
#
# Only rows where group == "treatment" will be processed.
#
# Expected peakcalling_params.tsv:
#   mark        format   genome   keep_dup   extra_params
#   H3K27ac     BAMPE    hs       1          "--call-summits -B --SPMR"
#   H3K4me3     BAMPE    hs       1          "--call-summits -B --SPMR"
#   H3K27me3    BAMPE    hs       1          "--broad --broad-cutoff 0.1 -B --SPMR"
#
# Outputs:
#   <out_dir>/<mark>/<sample_id>_peaks.*  (narrowPeak or broadPeak)
#   plus MACS2 logs and bedGraph files when applicable
#
# --------------------------------------------------------------------

SAMPLES_TSV="$1"
PARAMS_TSV="$2"
BAM_DIR="$3"
OUT_DIR="$4"

mkdir -p "${OUT_DIR}"

# --------------------------------------------------------------------
# Load per-mark MACS2 parameters into associative arrays
# --------------------------------------------------------------------
declare -A FORMAT
declare -A GENOME
declare -A KEEP_DUP
declare -A EXTRA

while IFS=$'\t' read -r mark format genome keep_dup extra_params; do
    [[ "$mark" == "mark" ]] && continue
    FORMAT["$mark"]="$format"
    GENOME["$mark"]="$genome"
    KEEP_DUP["$mark"]="$keep_dup"
    EXTRA["$mark"]="$extra_params"
done < "$PARAMS_TSV"

# --------------------------------------------------------------------
# Iterate through all treatment samples and call MACS2
# --------------------------------------------------------------------
while IFS=$'\t' read -r sample_id fastq1 fastq2 mark group control_sample; do
    [[ "$sample_id" == "sample_id" ]] && continue

    # Only call peaks for treatment samples
    if [[ "$group" != "treatment" ]]; then
        continue
    fi

    treat_bam="${BAM_DIR}/${sample_id}.bam"
    ctrl_bam="${BAM_DIR}/${control_sample}.bam"

    if [[ ! -f "$treat_bam" ]]; then
        echo "[WARN] Treatment BAM not found: $treat_bam"
        continue
    fi

    if [[ ! -f "$ctrl_bam" ]]; then
        echo "[WARN] Control BAM not found: $ctrl_bam"
        continue
    fi

    # Resolve MACS2 parameters for this histone mark
    macs_format="${FORMAT[$mark]:-BAMPE}"
    macs_genome="${GENOME[$mark]:-hs}"
    macs_keepdup="${KEEP_DUP[$mark]:-1}"
    macs_extra="${EXTRA[$mark]:-}"

    out_dir_mark="${OUT_DIR}/${mark}"
    mkdir -p "$out_dir_mark"

    echo "------------------------------------------------------------"
    echo "[INFO] Calling MACS2 peaks:"
    echo "       sample       = ${sample_id}"
    echo "       mark         = ${mark}"
    echo "       treatment    = ${treat_bam}"
    echo "       control      = ${ctrl_bam}"
    echo "       format       = ${macs_format}"
    echo "       genome       = ${macs_genome}"
    echo "       keep-dup     = ${macs_keepdup}"
    echo "       extra        = ${macs_extra}"
    echo "       output dir   = ${out_dir_mark}"
    echo "------------------------------------------------------------"

    macs2 callpeak \
        -t "$treat_bam" \
        -c "$ctrl_bam" \
        -f "$macs_format" \
        -g "$macs_genome" \
        --keep-dup "$macs_keepdup" \
        $macs_extra \
        -n "$sample_id" \
        --outdir "$out_dir_mark"

    echo "[INFO] Finished MACS2 for sample: $sample_id"
    echo
done < "$SAMPLES_TSV"
