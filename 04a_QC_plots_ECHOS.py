#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
04a_QC_plots_ECHOS_H3K4me3.py

QC and method-comparison plots for the ECHOS H3K4me3 benchmarking dataset.

This script is designed to be:
  - Reproducible and publication-ready
  - General enough to be adapted to other projects / marks
  - Paired with the command-line pipeline:
      01_trim_align_libcomplex.sh
      02_macs2_callpeaks.sh
      03_make_bigwig_and_SES.sh
      04_peakset_integration.sh

Inputs (typical):
  1) multiBigwigSummary output:
       - signal_matrix.npz  (for correlation heatmap)
  2) FRiP summary table:
       - FRiP_summary.tsv
  3) Peak BEDs for ECHOS vs ChIP:
       - ECHOS merged peaks
       - ChIP merged peaks

Outputs:
  - Correlation heatmap of per-sample signals
  - Boxplot of FRiP values by method
  - Area-proportional Venn diagram for ECHOS vs ChIP

All figures are saved as SVG / PDF / high-DPI PNG.
"""

from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


# -------------------------------------------------------------------------
# Global configuration
# -------------------------------------------------------------------------

# Output directory for all figures
OUTDIR = Path("pubfigs_ECHOS_H3K4me3")
OUTDIR.mkdir(parents=True, exist_ok=True)

# Files produced by 04_peakset_integration.sh (or equivalent)
FRIP_FILE = Path("results/peak_qc/FRiP_summary.tsv")
NPZ_CORR = Path("results/peak_qc/signal_matrix.npz")

# Peak BEDs for ECHOS vs ChIP Venn diagram
BED_ECHOS = Path("results/peak_qc/PSS_K4me3_2replicates_intersect_merge.bed")
BED_CHIP = Path("results/peak_qc/ChIP_K4me3_2replicates_intersect_merge.bed")

# Matplotlib defaults
matplotlib.rcParams.update(
    {
        "font.size": 10,
        "pdf.fonttype": 42,  # embed fonts in a vector-friendly way
        "svg.fonttype": "none",
    }
)

# Colors for methods (can be adapted)
METHOD_COLORS = {
    "ChIP": "#1f77b4",
    "ECHOS": "#ff7f0e",
    "CUTTag": "#2ca02c",
    "IgG": "#7f7f7f",
}

# Colormap for correlation heatmap
CMAP_CORR = "OrRd"

# Raster DPI for PNG exports
DPI_RASTER = 600


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

def save_figure(fig: matplotlib.figure.Figure, outbase: Path) -> None:
    """
    Save a Matplotlib figure to SVG, PDF, and PNG.
    outbase: path WITHOUT extension; extensions are added automatically.
    """
    outbase = Path(outbase)
    fig.savefig(outbase.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".png"), dpi=DPI_RASTER, bbox_inches="tight")
    print(f"Saved: {outbase}.svg / .pdf / .png")


# -------------------------------------------------------------------------
# 1) Load deepTools multiBigwigSummary outputs and plot correlation
# -------------------------------------------------------------------------

def load_multibigwigsummary_npz(npz_path: Path):
    """
    Load deepTools multiBigwigSummary output (.npz).

    Expected keys:
      - 'matrix': 2D array (regions x samples)
      - 'labels': 1D array of sample labels

    Returns:
      matrix: numpy array of shape (regions, samples)
      labels: list of sample labels (str)
    """
    arr = np.load(npz_path, allow_pickle=True)
    if "matrix" not in arr:
        # Fallback: try to infer matrix key
        matrix = None
        for k in arr.files:
            if isinstance(arr[k], np.ndarray) and arr[k].ndim == 2:
                matrix = arr[k]
                break
        if matrix is None:
            raise ValueError(
                f"No 2D 'matrix' found in {npz_path}. "
                f"Available keys: {list(arr.files)}"
            )
    else:
        matrix = arr["matrix"]

    if "labels" in arr:
        labels = [str(x) for x in arr["labels"]]
    else:
        labels = [f"sample_{i+1}" for i in range(matrix.shape[1])]

    return matrix, labels


def _clean_sample_label(raw: str) -> str:
    """
    Clean and standardize sample labels for plotting.

    This is project-specific and can be adapted to other datasets.
    For ECHOS H3K4me3, we:
      - strip '.bw' and '_treat_pileup'
      - map prefixes like 'ChIP_H3K4me3_' -> 'ChIP_'
      - map 'K4me3_' -> 'ECHOS_'
      - make 'rep1' -> 'rep 1' for better readability
    """
    x = raw
    x = x.replace(".bw", "").replace("_treat_pileup", "")
    x = re.sub(r"^ChIP_H3K4me3_", "ChIP_", x)
    x = re.sub(r"^ChIP_H3K4me3", "ChIP", x)
    x = re.sub(r"^K4me3_", "ECHOS_", x)
    x = x.replace("rep", "rep ")
    return x


def _infer_method(label_clean: str) -> str:
    """
    Infer method category from cleaned label.
    This is used to order samples and color-code them.
    """
    lc = label_clean.lower()
    if lc.startswith("chip"):
        return "ChIP"
    if lc.startswith("echos"):
        return "ECHOS"
    if lc.startswith("cuttag"):
        return "CUTTag"
    if "igg" in lc:
        return "IgG"
    return "Other"


def plot_correlation_heatmap(
    npz_path: Path = NPZ_CORR,
    out_prefix: str = "Fig_CorrHeatmap_ECHOS_H3K4me3",
) -> None:
    """
    Plot Pearson correlation heatmap from deepTools multiBigwigSummary output.

    Steps:
      1) Load matrix (regions x samples)
      2) Remove rows with NaN or all-zero signal
      3) Compute Pearson correlation between sample vectors
      4) Clean sample labels and order by method: ChIP -> ECHOS -> CUTTag -> IgG -> Other
      5) Plot heatmap with annotations and simple group separators
    """
    if not npz_path.is_file():
        print(f"[WARN] NPZ file not found, skip correlation heatmap: {npz_path}")
        return

    matrix, labels = load_multibigwigsummary_npz(npz_path)
    X = np.array(matrix, dtype=float)  # shape: (regions, samples)

    # Filter out rows with NaN or trivial signal
    row_mask = np.isfinite(X).all(axis=1) & (X.sum(axis=1) != 0)
    Xf = X[row_mask]

    # Compute correlation matrix (samples x samples)
    corr = np.corrcoef(Xf.T)

    # Clean labels and infer method categories
    labels_clean = [_clean_sample_label(l) for l in labels]
    methods = [_infer_method(lc) for lc in labels_clean]

    # Define preferred method order
    order_chip = [i for i, m in enumerate(methods) if m == "ChIP"]
    order_echos = [i for i, m in enumerate(methods) if m == "ECHOS"]
    order_cuttag = [i for i, m in enumerate(methods) if m == "CUTTag"]
    order_igg = [i for i, m in enumerate(methods) if m == "IgG"]
    order_other = [i for i, m in enumerate(methods) if m not in ("ChIP", "ECHOS", "CUTTag", "IgG")]

    order = order_chip + order_echos + order_cuttag + order_igg + order_other
    labels_ord = [labels_clean[i] for i in order]
    C = corr[np.ix_(order, order)]

    # Plot
    fig = plt.figure(figsize=(6.5, 5.8))
    ax = fig.add_subplot(111)

    im = ax.imshow(C, vmin=0, vmax=1, cmap=CMAP_CORR, aspect="equal")

    ax.set_xticks(range(len(labels_ord)))
    ax.set_yticks(range(len(labels_ord)))
    ax.set_xticklabels(labels_ord, rotation=40, ha="right")
    ax.set_yticklabels(labels_ord)

    # Annotate each cell with r (optional, but nice for publication-level QC)
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            val = C[i, j]
            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                fontsize=8,
                color=("white" if val >= 0.7 else "black"),
            )

    # Minor grid lines
    ax.set_xticks(np.arange(-0.5, len(labels_ord), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(labels_ord), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Group separators between ChIP / ECHOS / CUTTag / IgG / Other
    cuts = []
    offset = 0
    for group in [order_chip, order_echos, order_cuttag, order_igg, order_other]:
        if group:
            offset += len(group)
            cuts.append(offset - 0.5)
    # Don't draw the last cut at the end of matrix
    for c in cuts[:-1]:
        ax.axhline(c, color="black", linewidth=1.0)
        ax.axvline(c, color="black", linewidth=1.0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson r")
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

    ax.set_title("Correlation heatmap — ECHOS H3K4me3 benchmarking")
    fig.tight_layout()

    save_figure(fig, OUTDIR / out_prefix)
    plt.close(fig)
    print("Correlation heatmap done.")


# -------------------------------------------------------------------------
# 2) FRiP boxplot by method
# -------------------------------------------------------------------------

def _parse_method_from_sample(sample: str) -> str:
    """
    Infer method category from sample ID in FRiP_summary.tsv.
    This is project-specific; adjust patterns as needed.

    Examples:
      - CUTTag_... -> CUTTag
      - PSS_... or ECHOS_... -> ECHOS
      - ChIP_... -> ChIP
      - IgG_... -> IgG
    """
    s = sample.lower()
    if s.startswith("chip"):
        return "ChIP"
    if s.startswith("echos") or "pss" in s:
        return "ECHOS"
    if "cutntag" in s or s.startswith("cuttag"):
        return "CUTTag"
    if "igg" in s:
        return "IgG"
    return "Other"


def plot_frip_boxplot(
    frip_path: Path = FRIP_FILE,
    out_prefix: str = "Fig_FRiP_boxplot_ECHOS_H3K4me3",
) -> None:
    """
    Plot FRiP values (Fraction of Reads in Peaks) grouped by method.

    Input:
      FRiP_summary.tsv from 04_peakset_integration.sh, format:
        sample  total_frag_or_read   in_peaks    FRiP
    """
    if not frip_path.is_file():
        print(f"[WARN] FRiP file not found, skip FRiP boxplot: {frip_path}")
        return

    df = pd.read_csv(frip_path, sep="\t")
    if "FRiP" not in df.columns:
        print(f"[WARN] Column 'FRiP' not found in {frip_path}, skip FRiP plot.")
        return

    df["method"] = df["sample"].map(_parse_method_from_sample)

    # Sort methods in a stable, logical order
    method_order = ["ChIP", "ECHOS", "CUTTag", "IgG", "Other"]
    df["method"] = pd.Categorical(df["method"], categories=method_order, ordered=True)
    df = df.sort_values("method")

    fig, ax = plt.subplots(figsize=(4.5, 4.0))

    # Custom boxplot by method
    methods_present = [m for m in method_order if m in df["method"].unique()]
    positions = np.arange(len(methods_present))

    data_by_method = [df.loc[df["method"] == m, "FRiP"].values for m in methods_present]

    box = ax.boxplot(
        data_by_method,
        positions=positions,
        widths=0.6,
        patch_artist=True,
        showfliers=False,
    )

    # Color and styling
    for patch, m in zip(box["boxes"], methods_present):
        patch.set_facecolor(METHOD_COLORS.get(m, "#bbbbbb"))
        patch.set_edgecolor("black")
        patch.set_linewidth(1.0)
    for element in ["whiskers", "caps", "medians"]:
        for line in box[element]:
            line.set_color("black")
            line.set_linewidth(1.0)

    # Overlay individual points (jitter)
    for i, m in enumerate(methods_present):
        y = df.loc[df["method"] == m, "FRiP"].values
        x = np.random.normal(loc=positions[i], scale=0.06, size=len(y))
        ax.scatter(x, y, s=20, alpha=0.8, edgecolors="black", linewidths=0.3,
                   facecolors="white")

    ax.set_xticks(positions)
    ax.set_xticklabels(methods_present, rotation=0)
    ax.set_ylabel("FRiP (fraction of reads in peaks)")
    ax.set_title("FRiP by method — ECHOS H3K4me3")

    ax.set_ylim(0, min(1.0, max(0.05, df["FRiP"].max() * 1.2)))
    fig.tight_layout()

    save_figure(fig, OUTDIR / out_prefix)
    plt.close(fig)
    print("FRiP boxplot done.")


# -------------------------------------------------------------------------
# 3) Area-proportional Venn diagram for ECHOS vs ChIP peaks
# -------------------------------------------------------------------------

def read_bed(bed_path: Path):
    """
    Read a BED file into a dict:
      { chrom: [(start, end), ...], ... }
    Only uses first three columns (chrom, start, end).
    """
    bed_path = Path(bed_path)
    data = {}
    with bed_path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            xs = line.rstrip().split("\t")
            if len(xs) < 3:
                continue
            chrom, s, e = xs[0], int(xs[1]), int(xs[2])
            if s > e:
                s, e = e, s
            data.setdefault(chrom, []).append((s, e))
    for chrom in data:
        data[chrom].sort()
    return data


def merge_union_with_labels(bedA: dict, bedB: dict):
    """
    Union-based clustering of intervals from A/B, labeling each union cluster as:
      - A-only
      - B-only
      - Both

    Returns:
      a_only, b_only, both  (counts of union clusters in each category)
    """
    a_only = b_only = both = 0
    chroms = sorted(set(bedA.keys()) | set(bedB.keys()))

    for chrom in chroms:
        ivals = []
        for s, e in bedA.get(chrom, []):
            ivals.append((s, e, "A"))
        for s, e in bedB.get(chrom, []):
            ivals.append((s, e, "B"))
        if not ivals:
            continue
        ivals.sort()

        cur_s = cur_e = None
        hasA = hasB = False

        for s, e, lab in ivals:
            if cur_s is None:
                # start a new cluster
                cur_s, cur_e = s, e
                hasA = (lab == "A")
                hasB = (lab == "B")
                continue
            if s <= cur_e:
                # overlap with current cluster, extend
                if e > cur_e:
                    cur_e = e
                if lab == "A":
                    hasA = True
                else:
                    hasB = True
            else:
                # close current cluster
                if hasA and hasB:
                    both += 1
                elif hasA:
                    a_only += 1
                elif hasB:
                    b_only += 1
                # start new
                cur_s, cur_e = s, e
                hasA = (lab == "A")
                hasB = (lab == "B")

        # close last cluster
        if cur_s is not None:
            if hasA and hasB:
                both += 1
            elif hasA:
                a_only += 1
            elif hasB:
                b_only += 1

    return a_only, b_only, both


def _circle_overlap_area(r1: float, r2: float, d: float) -> float:
    """
    Compute area of overlap of two circles with radii r1, r2 and center distance d.
    """
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return np.pi * min(r1, r2) ** 2

    alpha = 2 * np.arccos((d ** 2 + r1 ** 2 - r2 ** 2) / (2 * d * r1))
    beta = 2 * np.arccos((d ** 2 + r2 ** 2 - r1 ** 2) / (2 * d * r2))
    a1 = 0.5 * r1 ** 2 * (alpha - np.sin(alpha))
    a2 = 0.5 * r2 ** 2 * (beta - np.sin(beta))
    return a1 + a2


def _solve_distance_for_overlap(r1: float, r2: float, target_area: float) -> float:
    """
    Numerical bisection to find distance d that yields desired overlap area.
    """
    lo, hi = 0.0, r1 + r2
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        area = _circle_overlap_area(r1, r2, mid)
        if area > target_area:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def plot_venn_two_sets(
    bed_a: Path = BED_ECHOS,
    bed_b: Path = BED_CHIP,
    label_a: str = "ECHOS",
    label_b: str = "ChIP",
    out_prefix: str = "Fig_Venn_ECHOS_vs_ChIP",
) -> None:
    """
    Area-proportional two-set Venn diagram for peak overlaps.

    The union clusters are counted as:
      - A-only
      - B-only
      - Both

    The circle areas are proportional to set sizes, and distance between
    circle centers is solved numerically to approximate the intersection area.
    """
    if (not bed_a.is_file()) or (not bed_b.is_file()):
        print(f"[WARN] BED files not found, skip Venn: {bed_a}, {bed_b}")
        return

    bedA = read_bed(bed_a)
    bedB = read_bed(bed_b)

    a_only, b_only, both = merge_union_with_labels(bedA, bedB)

    A = a_only + both
    B = b_only + both
    I = both
    U = a_only + b_only + both

    print(f"{label_a} total clusters: {A}")
    print(f"{label_b} total clusters: {B}")
    print(f"Intersection clusters:   {I}  |  Union: {U}")
    print(f"{label_a}-only: {a_only}   {label_b}-only: {b_only}   Both: {both}")

    # Radii proportional to sqrt(set size); rescale for nicer plot size
    rA0 = np.sqrt(A / np.pi) if A > 0 else 0.0
    rB0 = np.sqrt(B / np.pi) if B > 0 else 0.0
    max_r0 = max(rA0, rB0, 1e-9)
    target_max_r = 1.5
    k = (target_max_r / max_r0) ** 2
    r1 = np.sqrt(k * A / np.pi)
    r2 = np.sqrt(k * B / np.pi)
    target_overlap = k * I

    d = _solve_distance_for_overlap(r1, r2, target_overlap)

    fig, ax = plt.subplots(figsize=(6, 5))
    x1, x2 = -d / 2.0, d / 2.0

    # Colors: use a light and dark shade (Oranges palette-like)
    circle1 = plt.Circle(
        (x1, 0),
        r1,
        edgecolor="black",
        facecolor="#fdd0a2",  # light orange
        linewidth=1.5,
    )
    circle2 = plt.Circle(
        (x2, 0),
        r2,
        edgecolor="black",
        facecolor="#e6550d",  # dark orange
        alpha=0.55,
        linewidth=1.5,
    )
    ax.add_patch(circle1)
    ax.add_patch(circle2)

    # Percentages
    pctA = (I / A * 100.0) if A else 0.0
    pctB = (I / B * 100.0) if B else 0.0

    # Text annotations
    ax.text(
        x1 - r1 * 0.4,
        0,
        f"{label_a}-only\n{a_only}",
        ha="center",
        va="center",
    )
    ax.text(
        x2 + r2 * 0.4,
        0,
        f"{label_b}-only\n{b_only}",
        ha="center",
        va="center",
    )
    ax.text(
        0,
        0,
        f"Both\n{I}\n({pctA:.1f}% of {label_a}, {pctB:.1f}% of {label_b})",
        ha="center",
        va="center",
    )

    ax.set_aspect("equal")
    ax.set_xlim(min(x1 - r1, x2 - r2) - 0.3, max(x1 + r1, x2 + r2) + 0.3)
    ax.set_ylim(-max(r1, r2) - 0.3, max(r1, r2) + 0.8)
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_title(f"Area-proportional Venn — {label_a} vs {label_b}")
    ax.text(
        0,
        max(r1, r2) + 0.55,
        f"{label_a}={A}, {label_b}={B}, Intersection={I}, Union={U}",
        ha="center",
        va="center",
        fontsize=9,
    )

    fig.tight_layout()
    save_figure(fig, OUTDIR / out_prefix)
    plt.close(fig)
    print("Venn diagram done.")


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

def main():
    # 1) Correlation heatmap from multiBigwigSummary
    plot_correlation_heatmap(NPZ_CORR)

    # 2) FRiP boxplot by method
    plot_frip_boxplot(FRIP_FILE)

    # 3) Area-proportional Venn diagram (ECHOS vs ChIP)
    plot_venn_two_sets(
        bed_a=BED_ECHOS,
        bed_b=BED_CHIP,
        label_a="ECHOS",
        label_b="ChIP",
    )


if __name__ == "__main__":
    main()
