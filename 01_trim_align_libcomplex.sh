#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------------------------
# 01_trim_align_libcomplex.sh
#
# Universal preprocessing script for CUT&Tag / ECHOS datasets.
# Performs:
#   (1) Adapter trimming using BBDuk
#   (2) Alignment with Bowtie2
#   (3) Sorted BAM generation with samtools
#   (4) Library complexity estimation using Picard
#
# This script is intentionally marker-agnostic and works for any
# histone modification (H3K27ac, H3K4me3, H3K27me3, etc.).
#
# Usage:
#   bash scripts/01_trim_align_libcomplex.sh \
#       <fastq_dir> \
#       <output_dir> \
#       <bowtie2_index_prefix> \
#       <adapters_fasta> \
#       <picard_jar> \
#       <threads>
#
# Example:
#   bash scripts/01_trim_align_libcomplex.sh \
#       data/fastq \
#       results/align \
#       /indices/hg19/hg19 \
#       config/adapters.fa \
#       /opt/picard/picard.jar \
#       16
#
# Output structure:
#   <output_dir>/
#       bam/                 # Aligned and sorted BAM files
#       qc/                  # Library complexity metrics
#       *.trimmed_R*.fastq.gz  # Trimmed FASTQ files
#
# Requirements:
#   - bbduk.sh (BBMap)
#   - bowtie2
#   - samtools
#   - Picard (EstimateLibraryComplexity)
#   - Java 8+
#
# --------------------------------------------------------------------

FASTQ_DIR="$1"
OUT_DIR="$2"
BT2_INDEX="$3"
ADAPTERS="$4"
PICARD_JAR="$5"
THREADS="$6"

mkdir -p "${OUT_DIR}/bam" "${OUT_DIR}/qc"

# Iterate over all read pairs in the FASTQ directory
for fastq1 in "${FASTQ_DIR}"/*_R1_001.fastq.gz; do
    base=$(basename "${fastq1}" _R1_001.fastq.gz)
    fastq2="${FASTQ_DIR}/${base}_R2_001.fastq.gz"

    echo "------------------------------------------------------------------"
    echo "[INFO] Processing sample: ${base}"
    echo "------------------------------------------------------------------"

    # -----------------------------
    # 1) Adapter trimming (BBDuk)
    # -----------------------------
    trimmed1="${OUT_DIR}/${base}.trimmed_R1.fastq.gz"
    trimmed2="${OUT_DIR}/${base}.trimmed_R2.fastq.gz"

    echo "[INFO] Trimming adapters with BBDuk..."
    bbduk.sh -Xmx4g \
        in1="${fastq1}" \
        in2="${fastq2}" \
        out1="${trimmed1}" \
        out2="${trimmed2}" \
        ref="${ADAPTERS}" \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo

    # -----------------------------
    # 2) Alignment using Bowtie2
    # -----------------------------
    bam="${OUT_DIR}/bam/${base}.bam"

    echo "[INFO] Aligning with bowtie2..."
    bowtie2 -t -q \
        -N 1 -L 25 -X 2000 \
        -p "${THREADS}" \
        --local --no-mixed --no-discordant \
        -x "${BT2_INDEX}" \
        -1 "${trimmed1}" \
        -2 "${trimmed2}" \
        | samtools view -u \
        | samtools sort -o "${bam}"

    samtools index "${bam}"

    # --------------------------------------
    # 3) Library complexity (Picard)
    # --------------------------------------
    metrics="${OUT_DIR}/qc/${base}_lib_complexity.txt"

    echo "[INFO] Estimating library complexity (Picard)..."
    java -jar "${PICARD_JAR}" EstimateLibraryComplexity \
        I="${bam}" \
        O="${metrics}"

    echo "[INFO] Finished sample: ${base}"
    echo
done
