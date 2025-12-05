#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
07_micronuclei_chrY_hotspot_analysis.py

Micronuclei H3K27ac – chrY hotspot pipeline (generic version).

What this script does (high-level):
-----------------------------------
1) Configuration
   - Define project root, BAM locations and sample pairs
   - Create output directories

2) chrY QC (Part 4)
   - For each sample, compute per-chromosome coverage (reads / kb)
   - Compute chrY / median(autosomes) enrichment ratio
   - Save CSV + barplot (PDF)

3) Per-chromosome specificity (Part 5)
   - For micronuclei vs genomic reference:
     * Compute, for each chromosome, the proportion of reads in that chromosome
       relative to total mapped reads
     * Fold change = (MN proportion) / (genomic proportion)
   - Save CSV + multiple PDF plots (bar / scatter / histogram + summary)

4) SES-normalized chrY log2(IP/CTL) and hotspots
   - For each IP/control pair:
       bamCompare (deepTools) on chrY with SES, bin=5kb, smoothLength=50kb
   - multiBigwigSummary on chrY, 5kb bins for [MN1, MN2, TR, UT]
   - Compute:
       R_MN = mean(MN1, MN2)
       dR_TR = R_MN - TR
   - Call MN-specific domains:
       R_MN > 0.6 and dR_TR > 0.6
       Merge bins with max gap 5kb, keep domains >= 30kb
   - Call replicate-consistent domains:
       MN1 > 0.6, MN2 > 0.6, dR_TR > 0.6
       Same merging and length filtering
   - Plot SES log2 tracks with domain highlights

Dependencies:
-------------
Python: pysam, numpy, pandas, matplotlib, seaborn
CLI:   samtools, bedtools, deepTools (bamCoverage, bigwigCompare, multiBigwigSummary)
Modules: this script assumes an environment with Lmod/Environment Modules,
         but you can replace the `module load` lines with hard-coded paths.

NOTE:
-----
- Edit the CONFIG section below for your own project:
  * PROJECT_ROOT
  * MN_SAMPLES
  * WHOLE_SAMPLES
  * PAIRS_FOR_QC / PAIRS_FOR_SES
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import seaborn as sns


# ======================
# 0) CONFIGURATION
# ======================

CONFIG = {
    # Project root folder
    "PROJECT_ROOT": "/work/GCRB/s219796/PSS/MicroNuclei_202508",

    # BAM locations (edit to match your own structure)
    "MN_BAM_DIR": "analysis_T2T/1_aligned_bams",
    "WHOLE_BAM_DIR": "analysis_T2T_new_merged/1_aligned_bams",

    # Sample definitions
    # Micronuclei / nuclear H3K27ac & DNA (keys are arbitrary labels used below)
    "MN_SAMPLES": {
        "MN_K27ac_pss_1": "MN_K27ac_pss_1_T2T.bam",
        "MN_K27ac_pss_2": "MN_K27ac_pss_2_T2T.bam",
        "MN_whole_pss_1": "MN_whole_pss_1_T2T.bam",
        "MN_whole_pss_2": "MN_whole_pss_2_T2T.bam",
    },
    "WHOLE_MERGED_SAMPLES": {
        "K27ac_untreated_merged_UV1": "K27ac_untreated_UV1_merged_T2T.bam",
        "K27ac_treated_merged_UV1": "K27ac_treated_UV1_merged_T2T.bam",
        "whole_untreated_merged_UV1": "whole_untreated_UV1_merged_T2T.bam",
        "whole_untreated_merged_UV2": "whole_untreated_UV2_merged_T2T.bam",
        "whole_treated_merged_UV1": "whole_treated_UV1_merged_T2T.bam",
        "whole_treated_merged_UV2": "whole_treated_UV2_merged_T2T.bam",
    },

    # Pairs for QC (Part 4 / 5) and SES analysis:
    # IP samples (H3K27ac) and their nuclear controls
    "PAIRS_FOR_SES": [
        ("MN_K27ac_pss_1", "MN_whole_pss_1"),
        ("MN_K27ac_pss_2", "MN_whole_pss_2"),
        ("K27ac_untreated_merged_UV1", "whole_untreated_merged_UV1"),
        ("K27ac_treated_merged_UV1", "whole_treated_merged_UV1"),
    ],

    # Micronuclei group vs genomic reference used for per-chromosome specificity
    "MN_GROUP_FOR_SPECIFICITY": [
        "MN_K27ac_pss_1",
        "MN_K27ac_pss_2",
    ],
    "GENOMIC_GROUP_FOR_SPECIFICITY": [
        "K27ac_untreated_merged_UV1",
    ],

    # deepTools / bedtools / samtools threads
    "N_THREADS": 16,

    # Optional T2T blacklist (set to None if you don't want to use it)
    "BLACKLIST_BED": "T2T_blacklist.bed",
}


# ======================
# 1) HELPER FUNCTIONS
# ======================

def run_cmd_with_modules(cmd: str, extra_modules: List[str] = None) -> subprocess.CompletedProcess:
    """
    Run a shell command in a login-like Bash with Environment Modules.
    `extra_modules` is a list of module names to try `module load`.

    This is written in a generic way so it can be adapted to another cluster.
    """
    extra_modules = extra_modules or []
    module_init = (
        "source /etc/profile.d/modules.sh >/dev/null 2>&1 || "
        "source /usr/share/Modules/init/bash >/dev/null 2>&1 || true"
    )
    load_parts = []
    for m in extra_modules:
        load_parts.append(f"module load {m} >/dev/null 2>&1 || true")
    load_str = "; ".join(load_parts)
    full = f"{module_init}; {load_str}; {cmd}"
    print("$", cmd)
    proc = subprocess.run(
        ["bash", "-lc", full],
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr)
        raise RuntimeError(f"Command failed: {cmd}")
    return proc


def ensure_bai(bam: str) -> None:
    """Ensure BAM index exists (creates .bai if missing)."""
    bai = bam + ".bai"
    if os.path.exists(bai):
        return
    cmd = f"samtools index -@ {CONFIG['N_THREADS']} {bam}"
    run_cmd_with_modules(cmd, extra_modules=["samtools/1.10", "samtools"])


def get_all_bam_paths(project_root: Path) -> Dict[str, str]:
    """
    Build a dictionary {sample_key: absolute_bam_path} from CONFIG.
    """
    mn_dir = project_root / CONFIG["MN_BAM_DIR"]
    whole_dir = project_root / CONFIG["WHOLE_BAM_DIR"]

    paths = {}
    for k, rel in CONFIG["MN_SAMPLES"].items():
        paths[k] = str(mn_dir / rel)
    for k, rel in CONFIG["WHOLE_MERGED_SAMPLES"].items():
        paths[k] = str(whole_dir / rel)
    return paths


def load_chrom_lengths(bam_path: str) -> Dict[str, int]:
    """
    Read reference chromosome lengths from a BAM header.
    """
    lengths = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom in bam.header.references:
            lengths[chrom] = bam.get_reference_length(chrom)
    return lengths


def calculate_chromosome_coverage(bam_file_path: str,
                                  chrom_lengths: Dict[str, int]) -> Dict[str, float]:
    """
    Compute coverage per chromosome = read_count / (length_in_kb).

    This is used for chrY enrichment QC (Part 4).
    """
    coverage_data = {}
    try:
        with pysam.AlignmentFile(bam_file_path, "rb") as bamfile:
            for chrom_name, length in chrom_lengths.items():
                try:
                    read_count = bamfile.count(contig=chrom_name)
                    coverage = read_count / (length / 1000.0)
                    coverage_data[chrom_name] = coverage
                except ValueError:
                    # chromosome not present in this BAM
                    coverage_data[chrom_name] = 0.0
    except Exception as e:
        print(f"[WARN] Error reading BAM: {bam_file_path} ({e})")
        return {}
    return coverage_data


def calculate_raw_read_counts(bam_file_path: str,
                              chrom_lengths: Dict[str, int]) -> Dict[str, int]:
    """
    Compute raw read counts per chromosome.

    This is used for per-chromosome specificity (Part 5).
    """
    counts = {}
    try:
        with pysam.AlignmentFile(bam_file_path, "rb") as bamfile:
            for chrom_name in chrom_lengths.keys():
                try:
                    read_count = bamfile.count(contig=chrom_name)
                    counts[chrom_name] = read_count
                except ValueError:
                    counts[chrom_name] = 0
    except Exception as e:
        print(f"[WARN] Error reading BAM: {bam_file_path} ({e})")
        return {}
    return counts


def get_blacklist_arg(project_root: Path) -> str:
    """
    Return deepTools --blackListFileName argument if blacklist exists.
    """
    bl_rel = CONFIG.get("BLACKLIST_BED")
    if not bl_rel:
        return ""
    bl = project_root / bl_rel
    if bl.exists():
        return f"--blackListFileName {bl}"
    return ""


# ======================
# 2) PART 4 – chrY ENRICHMENT QC
# ======================

def run_chrY_enrichment_qc(project_root: Path,
                           bam_paths: Dict[str, str],
                           chrom_lengths: Dict[str, int],
                           out_root: Path) -> None:
    """
    Reproduce/simplify your Part 4:
    - For each sample, compute chrY coverage and median autosomal coverage.
    - Compute ratio = chrY / median_autosomes.
    - Save CSV + barplot.
    """
    print("--- Part 4: chrY enrichment QC ---")
    qc_dir = out_root / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # group structure can be flexibly changed here
    sample_groups = {
        "Micronuclei_H3K27ac": [
            ("MN_K27ac_Rep1", "MN_K27ac_pss_1"),
            ("MN_K27ac_Rep2", "MN_K27ac_pss_2"),
        ],
        "Micronuclei_DNA": [
            ("MN_whole_Rep1", "MN_whole_pss_1"),
            ("MN_whole_Rep2", "MN_whole_pss_2"),
        ],
        "Whole_Untreated_H3K27ac": [
            ("K27ac_Untreated_Merged", "K27ac_untreated_merged_UV1"),
        ],
        "Whole_Treated_H3K27ac": [
            ("K27ac_Treated_Merged", "K27ac_treated_merged_UV1"),
        ],
        "Whole_Untreated_DNA": [
            ("Whole_Untreated_UV1", "whole_untreated_merged_UV1"),
            ("Whole_Untreated_UV2", "whole_untreated_merged_UV2"),
        ],
        "Whole_Treated_DNA": [
            ("Whole_Treated_UV1", "whole_treated_merged_UV1"),
            ("Whole_Treated_UV2", "whole_treated_merged_UV2"),
        ],
    }

    rows = []
    for group_name, pairs in sample_groups.items():
        for sample_label, key in pairs:
            bam_path = bam_paths.get(key, "")
            if not os.path.exists(bam_path):
                print(f"[WARN] Missing BAM for {sample_label}: {bam_path}")
                continue
            cov = calculate_chromosome_coverage(bam_path, chrom_lengths)
            if not cov:
                continue
            chry_cov = cov.get("chrY", 0.0)
            autosome_covs = [cov.get(f"chr{i}", 0.0) for i in range(1, 23)]
            autosome_covs = [c for c in autosome_covs if c > 0.0]
            median_auto = float(np.median(autosome_covs)) if autosome_covs else 0.0
            ratio = chry_cov / (median_auto + 1e-6)
            rows.append({
                "Group": group_name,
                "Sample": sample_label,
                "chrY_Coverage": chry_cov,
                "Median_Autosome": median_auto,
                "Enrichment_Ratio": ratio,
            })

    df = pd.DataFrame(rows)
    out_csv = qc_dir / "y_chromosome_enrichment.csv"
    df.to_csv(out_csv, index=False)
    print("Saved:", out_csv)

    if not df.empty:
        plt.style.use("default")
        sns.set_palette("husl")
        plt.figure(figsize=(12, 6))
        sns.barplot(data=df, x="Sample", y="Enrichment_Ratio", hue="Group", dodge=False)
        plt.title("Y-Chromosome Enrichment", fontsize=16, weight="bold")
        plt.ylabel("chrY / median autosomes (coverage ratio)")
        plt.xlabel("Sample")
        plt.xticks(rotation=45, ha="right")
        plt.axhline(y=1.0, color="red", linestyle="--", alpha=0.5, label="No enrichment")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        out_pdf = qc_dir / "y_enrichment_barplot.pdf"
        plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
        plt.close()
        print("Saved:", out_pdf)


# ======================
# 3) PART 5 – PER-CHROMOSOME SPECIFICITY
# ======================

def run_per_chromosome_specificity(project_root: Path,
                                   bam_paths: Dict[str, str],
                                   chrom_lengths: Dict[str, int],
                                   out_root: Path) -> None:
    """
    Reproduce/simplify Part 5:
    - Compare micronuclei vs genomic H3K27ac per chromosome.
    - Compute proportion of reads for each chromosome in MN vs genomic.
    - Fold change = MN_prop / genomic_prop.
    - Save CSV + multiple summary plots (similar to your subplot layout).
    """
    print("--- Part 5: per-chromosome specificity ---")
    spec_dir = out_root / "part5_specificity"
    spec_dir.mkdir(parents=True, exist_ok=True)

    mn_keys = CONFIG["MN_GROUP_FOR_SPECIFICITY"]
    gen_keys = CONFIG["GENOMIC_GROUP_FOR_SPECIFICITY"]

    fold_results = []
    for chrom_name in sorted(chrom_lengths.keys()):
        if not chrom_name.startswith("chr"):
            continue

        mn_chr_reads, mn_tot_reads = [], []
        for k in mn_keys:
            bam_path = bam_paths.get(k, "")
            if not os.path.exists(bam_path):
                continue
            raw_counts = calculate_raw_read_counts(bam_path, chrom_lengths)
            mn_chr_reads.append(raw_counts.get(chrom_name, 0))
            mn_tot_reads.append(sum(raw_counts.values()))

        g_chr_reads, g_tot_reads = [], []
        for k in gen_keys:
            bam_path = bam_paths.get(k, "")
            if not os.path.exists(bam_path):
                continue
            raw_counts = calculate_raw_read_counts(bam_path, chrom_lengths)
            g_chr_reads.append(raw_counts.get(chrom_name, 0))
            g_tot_reads.append(sum(raw_counts.values()))

        if mn_chr_reads and g_chr_reads:
            avg_mn_chr = np.mean(mn_chr_reads)
            avg_mn_tot = np.mean(mn_tot_reads)
            avg_g_chr = np.mean(g_chr_reads)
            avg_g_tot = np.mean(g_tot_reads)
            mn_prop = avg_mn_chr / (avg_mn_tot + 1e-6)
            g_prop = avg_g_chr / (avg_g_tot + 1e-6)
            fc = mn_prop / (g_prop + 1e-6)
            fold_results.append({
                "Chromosome": chrom_name,
                "Micronuclei_Proportion": mn_prop,
                "Genomic_Proportion": g_prop,
                "Fold_Change": fc,
                "Log2FC": np.log2(fc),
            })

    fc_df = pd.DataFrame(fold_results)
    out_csv = spec_dir / "fold_change_specificity.csv"
    fc_df.to_csv(out_csv, index=False)
    print("Saved:", out_csv)

    if fc_df.empty:
        return

    # Plots (similar to your multi-panel figure)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle("Per-chromosome specificity (MN vs genomic)", fontsize=18, fontweight="bold")

    # 1) Fold-change barplot
    ax1 = axes[0, 0]
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    plot_data = fc_df.set_index("Chromosome").reindex(chrom_order).dropna()
    colors = ["red" if x > 1.5 else "blue" if x < 0.67 else "gray"
              for x in plot_data["Fold_Change"]]
    ax1.bar(range(len(plot_data)), plot_data["Fold_Change"], color=colors, alpha=0.7)
    ax1.set_xlabel("Chromosome")
    ax1.set_ylabel("Fold Change (MN / genomic)")
    ax1.set_title("Fold change by chromosome")
    ax1.set_xticks(range(len(plot_data)))
    ax1.set_xticklabels(plot_data.index, rotation=45, ha="right")
    ax1.axhline(y=1.0, color="black", linestyle="--", alpha=0.7)
    ax1.axhline(y=1.5, color="red", linestyle=":", alpha=0.5)
    ax1.axhline(y=0.67, color="blue", linestyle=":", alpha=0.5)
    ax1.grid(True, alpha=0.3)

    # 2) Proportion scatter
    ax2 = axes[0, 1]
    sc = ax2.scatter(fc_df["Genomic_Proportion"],
                     fc_df["Micronuclei_Proportion"],
                     c=fc_df["Fold_Change"],
                     cmap="RdYlBu_r",
                     s=80,
                     alpha=0.8)
    ax2.set_xlabel("Genomic proportion")
    ax2.set_ylabel("Micronuclei proportion")
    ax2.set_title("Proportion comparison")
    max_prop = max(fc_df["Genomic_Proportion"].max(),
                   fc_df["Micronuclei_Proportion"].max())
    ax2.plot([0, max_prop], [0, max_prop], "k--", alpha=0.5)
    ax2.grid(True, alpha=0.3)
    fig.colorbar(sc, ax=ax2, label="Fold change")

    # 3) Log2FC histogram
    ax3 = axes[1, 0]
    ax3.hist(fc_df["Log2FC"], bins=15, alpha=0.7, edgecolor="black")
    ax3.set_xlabel("Log2 Fold Change (MN / genomic)")
    ax3.set_ylabel("Count")
    ax3.set_title("Log2FC distribution")
    ax3.axvline(x=0, color="red", linestyle="--", alpha=0.7)
    ax3.axvline(x=1, color="orange", linestyle="--", alpha=0.7)
    ax3.axvline(x=-1, color="green", linestyle="--", alpha=0.7)
    ax3.grid(True, alpha=0.3)

    # 4) Summary text
    ax4 = axes[1, 1]
    ax4.axis("off")
    mean_fc = fc_df["Fold_Change"].mean()
    median_fc = fc_df["Fold_Change"].median()
    txt = (
        "Summary statistics\n\n"
        f"Mean FC:    {mean_fc:.3f}\n"
        f"Median FC:  {median_fc:.3f}\n"
        f"Min FC:     {fc_df['Fold_Change'].min():.3f}\n"
        f"Max FC:     {fc_df['Fold_Change'].max():.3f}\n"
        f"N chrom:    {len(fc_df)}\n"
    )
    ax4.text(0.05, 0.95, txt,
             transform=ax4.transAxes,
             va="top",
             family="monospace",
             fontsize=11,
             bbox=dict(boxstyle="round", facecolor="lightgray", alpha=0.85))

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_pdf = spec_dir / "fold_change_specificity_plots.pdf"
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", out_pdf)


# ======================
# 4) SES chrY log2(IP/CTL) & HOTSPOTS
# ======================

def build_ses_bigwigs_chrY(project_root: Path,
                           bam_paths: Dict[str, str],
                           out_root: Path) -> Dict[str, str]:
    """
    For each IP/control pair in CONFIG['PAIRS_FOR_SES']:

    Use deepTools bamCompare with SES normalization, restricted to chrY:
      - binSize 5000
      - smoothLength 50000
      - log2(IP/CTL)
    Returns dict {sample_key: bw_path}.
    """
    print("--- SES log2(IP/CTL) bigWigs on chrY (5kb / 50kb) ---")
    bw_dir = out_root / "chrY_only" / "bw"
    bw_dir.mkdir(parents=True, exist_ok=True)

    blacklist_arg = get_blacklist_arg(project_root)
    bw_paths = {}

    for ip_key, ct_key in CONFIG["PAIRS_FOR_SES"]:
        ip_bam = bam_paths.get(ip_key, "")
        ct_bam = bam_paths.get(ct_key, "")
        if not (os.path.exists(ip_bam) and os.path.exists(ct_bam)):
            print(f"[WARN] Skip SES (missing BAM): {ip_key}, {ct_key}")
            continue
        ensure_bai(ip_bam)
        ensure_bai(ct_bam)
        out_bw = bw_dir / f"{ip_key}.SES.5kb.smooth50kb.chrY.bw"
        cmd = (
            f"bamCompare -b1 {ip_bam} -b2 {ct_bam} "
            f"--operation log2 --pseudocount 1 "
            f"--region chrY --binSize 5000 --smoothLength 50000 "
            f"--scaleFactorsMethod SES --minMappingQuality 30 "
            f"-p {CONFIG['N_THREADS']} {blacklist_arg} "
            f"-o {out_bw}"
        )
        run_cmd_with_modules(cmd, extra_modules=["deeptools/3.5.0", "deeptools"])
        bw_paths[ip_key] = str(out_bw)
        print("Saved SES bw:", out_bw)

    return bw_paths


def ses_matrix_and_hotspots(project_root: Path,
                            bw_paths: Dict[str, str],
                            out_root: Path) -> None:
    """
    Use multiBigwigSummary on chrY (bin=5kb) with SES bigWigs for
    MN1, MN2, TR, UT, then:

    - Compute:
        R_MN = mean(MN1, MN2)
        dR_TR = R_MN - TR
    - Call MN-specific hotspots:
        R_MN > 0.6 and dR_TR > 0.6
      Merge bins with gap <= 5kb, keep domains >= 30kb.
    - Call replicate-consistent hotspots:
        MN1 > 0.6, MN2 > 0.6, dR_TR > 0.6
      Same merging / filtering.
    - Plot SES tracks with domain highlights.
    """
    print("--- SES matrix, domain-level hotspots, and plotting ---")

    # Keys used for matrix; you can rename to your own sample IDs if needed.
    MN1 = "MN_K27ac_pss_1"
    MN2 = "MN_K27ac_pss_2"
    TR = "K27ac_treated_merged_UV1"
    UT = "K27ac_untreated_merged_UV1"

    ses_dir = out_root / "chrY_only" / "bw" / "summary_ses"
    ses_dir.mkdir(parents=True, exist_ok=True)

    # Check that all four bigWigs are available
    needed = [MN1, MN2, TR, UT]
    if not all(k in bw_paths for k in needed):
        print("[WARN] Missing some SES bigWigs, skip SES hotspot analysis.")
        return

    bw_list = [bw_paths[MN1], bw_paths[MN2], bw_paths[TR], bw_paths[UT]]
    raw_tab = ses_dir / "chrY_log2_5kb50_smooth.tab"
    npz = ses_dir / "chrY_log2_5kb50_smooth.npz"

    cmd = (
        f"multiBigwigSummary bins --binSize 5000 --region chrY "
        f"-b {' '.join(bw_list)} "
        f"-p {CONFIG['N_THREADS']} "
        f"-out {npz} --outRawCounts {raw_tab}"
    )
    run_cmd_with_modules(cmd, extra_modules=["deeptools/3.5.0", "deeptools"])

    # Load matrix
    df = pd.read_csv(raw_tab, sep="\t", header=0, comment="#")
    # standardize columns: first three = coords, last four = samples
    df = df.rename(columns={
        df.columns[0]: "chromosome",
        df.columns[1]: "start",
        df.columns[2]: "end",
    })
    tails = list(df.columns)[-4:]
    df = df.rename(columns={
        tails[0]: "MN1",
        tails[1]: "MN2",
        tails[2]: "TR",
        tails[3]: "UT",
    })
    for c in ["start", "end", "MN1", "MN2", "TR", "UT"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["mid"] = (df["start"] + df["end"]) // 2

    # R_MN and dR_TR
    df["R_MN"] = df[["MN1", "MN2"]].mean(axis=1)
    df["dR_TR"] = df["R_MN"] - df["TR"]

    # Save full matrix
    R_matrix_csv = ses_dir / "chrY_log2_R_values_5kb50.csv"
    df.to_csv(R_matrix_csv, index=False)
    print("Saved:", R_matrix_csv)

    # --- MN-specific domains (threshold on R_MN and dR_TR) ---
    mask_mn = (df["R_MN"] > 0.6) & (df["dR_TR"] > 0.6)
    sel_mn = df.loc[mask_mn, ["chromosome", "start", "end", "R_MN"]]
    bins_bed = ses_dir / "MN_bins_5kb50.sel.bed"
    sel_mn.to_csv(bins_bed, sep="\t", header=False, index=False)

    merged_bed = ses_dir / "MN_domains_5kb50.merged.bed"
    cmd_merge = f"bedtools sort -i {bins_bed} | bedtools merge -d 5000 -i - -c 4 -o mean > {merged_bed}"
    run_cmd_with_modules(cmd_merge, extra_modules=["bedtools/2.30.0", "bedtools"])

    filt_bed = ses_dir / "MN_domains_5kb50.merged.min30kb.bed"
    with open(filt_bed, "w") as w, open(merged_bed) as r:
        for ln in r:
            ch, s, e, *rest = ln.strip().split("\t")
            s = int(s)
            e = int(e)
            if e - s >= 30000:
                w.write(ln)
    print("Saved MN-specific domains:", filt_bed)

    # --- Replicate-consistent domains (MN1 & MN2 both > 0.6, dR_TR > 0.6) ---
    mask_rep = (df["MN1"] > 0.6) & (df["MN2"] > 0.6) & (df["dR_TR"] > 0.6)
    sel_rep = df.loc[mask_rep, ["chromosome", "start", "end", "R_MN"]]
    bins_rep = ses_dir / "MN_bins_repConsistent_5kb50.sel.bed"
    sel_rep.to_csv(bins_rep, sep="\t", header=False, index=False)

    merged_rep = ses_dir / "MN_domains_repConsistent_5kb50.merged.bed"
    cmd_merge_rep = f"bedtools sort -i {bins_rep} | bedtools merge -d 5000 -i - -c 4 -o mean > {merged_rep}"
    run_cmd_with_modules(cmd_merge_rep, extra_modules=["bedtools/2.30.0", "bedtools"])

    filt_rep = ses_dir / "MN_domains_repConsistent_5kb50.merged.min30kb.bed"
    with open(filt_rep, "w") as w, open(merged_rep) as r:
        for ln in r:
            ch, s, e, *rest = ln.strip().split("\t")
            s = int(s)
            e = int(e)
            if e - s >= 30000:
                w.write(ln)
    print("Saved replicate-consistent domains:", filt_rep)

    # Quick summary
    def _summary(path: Path) -> Tuple[int, float]:
        n, L = 0, 0
        if not path.exists():
            return n, L
        with open(path) as r:
            for ln in r:
                ch, s, e, *rest = ln.strip().split("\t")
                s = int(s)
                e = int(e)
                n += 1
                L += (e - s)
        return n, L / 1e6

    n_dom, L_dom = _summary(filt_bed)
    n_rep, L_rep = _summary(filt_rep)
    print(f"MN-specific domains: {n_dom} (total length {L_dom:.2f} Mb)")
    print(f"Replicate-consistent domains: {n_rep} (total length {L_rep:.2f} Mb)")

    # --- Plot SES tracks with domain highlights ---
    print("--- Plot SES 5kb/50kb tracks with domain highlights ---")
    doms = []
    if filt_bed.exists():
        with open(filt_bed) as r:
            for ln in r:
                ch, s, e, *rest = ln.strip().split("\t")
                doms.append((int(s), int(e)))
    doms_rep = []
    if filt_rep.exists():
        with open(filt_rep) as r:
            for ln in r:
                ch, s, e, *rest = ln.strip().split("\t")
                doms_rep.append((int(s), int(e)))

    # Use df (same as matrix) for plotting
    x = df["mid"].to_numpy()
    plt.figure(figsize=(18, 4))
    plt.scatter(x, df["MN1"], s=6, label="MN1", color="#1f77b4", alpha=0.8)
    plt.scatter(x, df["MN2"], s=6, label="MN2", color="#ff7f0e", alpha=0.8)
    plt.scatter(x, df["TR"], s=6, label="TR", color="#2ca02c", alpha=0.8)
    plt.scatter(x, df["UT"], s=6, label="UT", color="#d62728", alpha=0.8)

    for s, e in doms:
        plt.axvspan(s, e, color="gold", alpha=0.15, lw=0)
    for s, e in doms_rep:
        plt.axvspan(s, e, color="red", alpha=0.10, lw=0)

    plt.title("chrY SES-smoothed log2(IP/CTL) (5kb bin, 50kb smooth)")
    plt.xlabel("chrY position (bp)")
    plt.ylabel("log2(IP/CTL)")
    plt.legend(ncol=4, fontsize=9)
    plt.grid(alpha=0.2, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    out_pdf1 = ses_dir / "overlay.SES.5kb50.track.pdf"
    plt.savefig(out_pdf1, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", out_pdf1)

    # MN-only with replicate-consistent domains
    plt.figure(figsize=(18, 4))
    plt.scatter(x, df["MN1"], s=6, label="MN1", color="#1f77b4", alpha=0.8)
    plt.scatter(x, df["MN2"], s=6, label="MN2", color="#ff7f0e", alpha=0.8)
    for s, e in doms_rep:
        plt.axvspan(s, e, color="red", alpha=0.12, lw=0)
    plt.title("chrY MN replicate-consistent domains (5kb bin, 50kb smooth)")
    plt.xlabel("chrY position (bp)")
    plt.ylabel("log2(IP/CTL)")
    plt.legend(fontsize=10)
    plt.grid(alpha=0.2, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    out_pdf2 = ses_dir / "MN_only.SES.5kb50.track.pdf"
    plt.savefig(out_pdf2, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", out_pdf2)


# ======================
# 5) MAIN
# ======================

def main():
    project_root = Path(CONFIG["PROJECT_ROOT"]).resolve()
    out_root = project_root / "analysis_T2T_new_merged" / "2_peak_calling"
    out_root.mkdir(parents=True, exist_ok=True)

    # Build BAM path dict and sanity check
    bam_paths = get_all_bam_paths(project_root)
    print("Project root:", project_root)
    print("Output root :", out_root)
    print("Samples & BAMs:")
    for k, v in sorted(bam_paths.items()):
        print(f"  {k}: {'OK' if os.path.exists(v) else 'MISSING'} - {v}")

    # Chromosome lengths (use one MN H3K27ac as reference)
    ref_key = CONFIG["MN_GROUP_FOR_SPECIFICITY"][0]
    ref_bam = bam_paths.get(ref_key, "")
    if not os.path.exists(ref_bam):
        raise FileNotFoundError(f"Reference BAM for chrom lengths not found: {ref_bam}")
    chrom_lengths = load_chrom_lengths(ref_bam)
    print(f"Loaded {len(chrom_lengths)} chromosomes from header of {ref_key}")

    # Part 4: chrY enrichment
    run_chrY_enrichment_qc(project_root, bam_paths, chrom_lengths, out_root)

    # Part 5: per-chromosome specificity
    run_per_chromosome_specificity(project_root, bam_paths, chrom_lengths, out_root)

    # SES chrY bigWigs and hotspots
    bw_paths = build_ses_bigwigs_chrY(project_root, bam_paths, out_root)
    ses_matrix_and_hotspots(project_root, bw_paths, out_root)


if __name__ == "__main__":
    main()
