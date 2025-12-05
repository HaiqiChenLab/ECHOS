#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
05_diff_peaks_prepare_counts.py

Generic script to prepare peak-by-sample count matrices for differential
peak analysis (e.g. with DESeq2 in 05a_diff_peaks_DESeq2.R).

What this script does:
  1) Takes per-sample peak files (e.g. MACS2 narrowPeak/broadPeak) and
     builds a union peak set (bed3) using bedtools merge.
  2) Uses bedtools multicov to count reads in each union peak for all samples.
  3) Writes:
       - counts_peak_matrix.csv (rows = peaks, columns = samples)
       - sample_metadata.csv (sample, group, bam, peaks)

Dependencies:
  - Python 3
  - bedtools (bedtools merge, bedtools multicov)
  - samtools (for BAM indexing if needed)

Usage:
  python 05_diff_peaks_prepare_counts.py

Edit the CONFIG section below to match your project.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import List, Dict

import pandas as pd

# ============================================================
# CONFIG (EDIT THIS BLOCK)
# ============================================================

# Base directory of the project (optional, can also use absolute paths below)
WORKDIR = Path("/path/to/project").resolve()

# Output directory for this contrast
OUTDIR = WORKDIR / "results" / "diff_peaks" / "example_contrast"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Distance (in bp) for merging peaks into union clusters
MERGE_DISTANCE = 1000

# List of samples for this contrast.
# Each entry must define:
#   - sample: sample ID used as column name in the count matrix
#   - bam:    path to the BAM file
#   - peaks:  path to the peak file (BED / narrowPeak / broadPeak; at least 3 columns: chr, start, end)
#   - group:  group label (e.g. "basal", "luminal", "young", "aged")
#
# Example below mimics basal vs luminal in cervix:
SAMPLES: List[Dict[str, str]] = [
    {
        "sample": "basal_rep1",
        "bam": str(WORKDIR / "bam" / "basal_rep1.bam"),
        "peaks": str(WORKDIR / "peaks" / "basal_rep1_peaks.narrowPeak"),
        "group": "basal",
    },
    {
        "sample": "basal_rep2",
        "bam": str(WORKDIR / "bam" / "basal_rep2.bam"),
        "peaks": str(WORKDIR / "peaks" / "basal_rep2_peaks.narrowPeak"),
        "group": "basal",
    },
    {
        "sample": "luminal_rep1",
        "bam": str(WORKDIR / "bam" / "luminal_rep1.bam"),
        "peaks": str(WORKDIR / "peaks" / "luminal_rep1_peaks.narrowPeak"),
        "group": "luminal",
    },
    {
        "sample": "luminal_rep3",
        "bam": str(WORKDIR / "bam" / "luminal_rep3.bam"),
        "peaks": str(WORKDIR / "peaks" / "luminal_rep3_peaks.narrowPeak"),
        "group": "luminal",
    },
]

# Filenames for outputs (relative to OUTDIR)
UNION_BED = OUTDIR / "union_peaks.bed"
COUNTS_CSV = OUTDIR / "counts_peak_matrix.csv"
META_CSV = OUTDIR / "sample_metadata.csv"

# ============================================================
# END OF CONFIG
# ============================================================


def say(*args):
    """Small logger."""
    print("[05_prepare_counts]", *args)


def ensure_tool(name: str):
    """Check that an external tool (bedtools, samtools) is available."""
    path = subprocess.getoutput(f"bash -lc 'which {name} || true'")
    if not path:
        raise RuntimeError(f"Required tool '{name}' not found in PATH.")
    return path


def build_union_peaks(samples: List[Dict[str, str]], union_bed: Path, merge_dist: int = 1000):
    """
    Build a union peak set from sample-specific peak files.

    Steps:
      1) Concatenate all peak files (chr, start, end) → tmp
      2) Sort + merge using bedtools → union_bed
    """
    if union_bed.exists():
        say("Union BED already exists, will overwrite:", union_bed)

    tmp_cat = union_bed.with_suffix(".tmp_cat.bed")
    tmp_sorted = union_bed.with_suffix(".tmp_sorted.bed")

    say("Collecting peak files and writing raw concatenated BED...")
    with tmp_cat.open("w") as w:
        n_lines = 0
        for s in samples:
            peak_path = Path(s["peaks"])
            if not peak_path.exists():
                raise FileNotFoundError(f"Peak file not found: {peak_path}")
            with peak_path.open() as r:
                for line in r:
                    if not line.strip() or line.startswith("#"):
                        continue
                    fs = line.rstrip().split("\t")
                    if len(fs) < 3:
                        continue
                    w.write("\t".join(fs[:3]) + "\n")
                    n_lines += 1
        say(f"  wrote {n_lines} lines from {len(samples)} peak files.")

    if n_lines == 0:
        raise RuntimeError("No valid peak lines found across all peak files.")

    # Ensure bedtools is available
    bedtools = ensure_tool("bedtools")

    say("Sorting concatenated peaks...")
    cmd_sort = f"{bedtools} sort -i {tmp_cat} > {tmp_sorted}"
    say("  ", cmd_sort)
    subprocess.run(["bash", "-lc", cmd_sort], check=True)

    say(f"Merging peaks with distance {merge_dist} bp...")
    cmd_merge = f"{bedtools} merge -i {tmp_sorted} -d {merge_dist} > {union_bed}"
    say("  ", cmd_merge)
    subprocess.run(["bash", "-lc", cmd_merge], check=True)

    tmp_cat.unlink(missing_ok=True)
    tmp_sorted.unlink(missing_ok=True)

    n_union = sum(1 for _ in union_bed.open())
    say(f"Union peaks written to {union_bed} ({n_union} intervals).")


def ensure_bam_indices(samples: List[Dict[str, str]]):
    """
    Ensure all BAM files have .bai indices (using samtools index if needed).
    """
    samtools = ensure_tool("samtools")

    for s in samples:
        bam = Path(s["bam"])
        if not bam.exists():
            raise FileNotFoundError(f"BAM not found: {bam}")
        bai = bam.with_suffix(bam.suffix + ".bai")
        if bai.exists():
            say("BAM index exists:", bai.name)
            continue
        cmd_idx = f"{samtools} index -@ 4 {bam}"
        say("Indexing BAM:", cmd_idx)
        subprocess.run(["bash", "-lc", cmd_idx], check=True)


def count_with_multicov(
    samples: List[Dict[str, str]],
    union_bed: Path,
    out_csv: Path,
):
    """
    Use bedtools multicov to count reads in each union peak for all samples.

    Output:
      - out_csv: CSV with rows = peaks (chr_start_end), columns = sample IDs.
    """
    bedtools = ensure_tool("bedtools")

    if not union_bed.exists():
        raise FileNotFoundError(f"Union BED does not exist: {union_bed}")

    # Ensure BAM indices
    ensure_bam_indices(samples)

    bams = [str(Path(s["bam"])) for s in samples]
    tmp_multicov = union_bed.with_suffix(".multicov.txt")

    say("Running bedtools multicov...")
    cmd_mc = f"{bedtools} multicov -bams {' '.join(bams)} -bed {union_bed} > {tmp_multicov}"
    say("  ", cmd_mc)
    subprocess.run(["bash", "-lc", cmd_mc], check=True)

    say("Parsing multicov output...")
    df = pd.read_csv(tmp_multicov, sep="\t", header=None)

    # Expect: chr, start, end, counts for each BAM
    min_cols = 3 + len(bams)
    if df.shape[1] < min_cols:
        raise RuntimeError(
            f"multicov output has fewer columns ({df.shape[1]}) than expected ({min_cols})."
        )

    chr_col = df.iloc[:, 0].astype(str)
    start_col = df.iloc[:, 1].astype(int).astype(str)
    end_col = df.iloc[:, 2].astype(int).astype(str)
    peak_ids = chr_col + "_" + start_col + "_" + end_col

    counts = df.iloc[:, 3 : 3 + len(bams)].copy()
    counts.index = peak_ids
    counts.columns = [s["sample"] for s in samples]

    # Save as CSV
    counts.to_csv(out_csv)
    tmp_multicov.unlink(missing_ok=True)

    say(f"Count matrix written to {out_csv} (peaks x samples = {counts.shape[0]} x {counts.shape[1]}).")


def write_sample_metadata(samples: List[Dict[str, str]], out_csv: Path):
    """
    Write sample_metadata.csv with columns:
      sample, group, bam, peaks
    """
    rows = []
    for s in samples:
        rows.append(
            {
                "sample": s["sample"],
                "group": s["group"],
                "bam": s["bam"],
                "peaks": s["peaks"],
            }
        )
    meta = pd.DataFrame(rows)
    meta.to_csv(out_csv, index=False)
    say(f"Sample metadata written to {out_csv}.")


def main():
    say("Starting 05_diff_peaks_prepare_counts.py")
    say("WORKDIR:", WORKDIR)
    say("OUTDIR :", OUTDIR)
    say("Number of samples:", len(SAMPLES))

    # 1) Build union peak set
    build_union_peaks(SAMPLES, UNION_BED, MERGE_DISTANCE)

    # 2) Count reads in union peaks for each sample
    count_with_multicov(SAMPLES, UNION_BED, COUNTS_CSV)

    # 3) Write sample metadata
    write_sample_metadata(SAMPLES, META_CSV)

    say("DONE. You can now run 05a_diff_peaks_DESeq2.R using:")
    say("  counts_file =", COUNTS_CSV)
    say("  meta_file   =", META_CSV)


if __name__ == "__main__":
    main()
