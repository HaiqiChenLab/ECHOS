#!/usr/bin/env python
"""
06_rna_ECHOS_integration.py

Generic RNA × ECHOS integration (e.g. RNA-seq vs promoter H3K4me3).

This script:

1) Ensures the union peak set has a gene annotation table
   (e.g. HOMER annotatePeaks output).
2) Aggregates chromatin (ECHOS / ChIP / CUT&Tag) DE results to gene level
   (best peak per gene; and, when available, promoter-only aggregation).
3) Loads RNA DESeq2 results and harmonizes gene identifiers
   (Ensembl ID preferred; gene symbol as optional fallback).
4) Joins RNA and chromatin on gene IDs.
5) Classifies genes into quadrants based on RNA and chromatin log2FC:
   - both up (Quadrant I)
   - both down (Quadrant III)
   - opposite directions
   - non-significant
6) Produces a scatter plot and exports gene lists and summary tables.

The logic is generic and can be reused for any RNA–chromatin contrast,
as long as the input formats follow the assumptions described below.
"""

from pathlib import Path
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

# ---------------------------------------------------------
# 0) Optional: Spearman correlation (scipy)
# ---------------------------------------------------------
try:
    from scipy.stats import spearmanr
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# =========================================================
# CONFIG (EDIT THIS BLOCK)
# =========================================================

# Project root (only used if you want to build paths relative to it)
WORKDIR = Path("/path/to/project_root")

# ---------------- Chromatin / ECHOS inputs ---------------

# Chromatin DE results (DESeq2 output on union peaks)
# Prefer unshrunken DESeq2 result if available; otherwise use the shrunken one.
CHIP_DE_DIR = WORKDIR / "results/chromatin_example/deseq2_results"
CHIP_DE_UNSHRUNKEN = CHIP_DE_DIR / "deseq2_peaks_contrast_unshrunken.csv"
CHIP_DE_SHRUNKEN   = CHIP_DE_DIR / "deseq2_peaks_contrast_shrunken.csv"

# Union peaks (BED3 format) used for DESeq2 on peaks
CHIP_UNION_BED = WORKDIR / "results/chromatin_example/union_peaks.bed"

# Peak annotation table (TSV) linking peak_id to genome + gene info.
# If this file does NOT exist and you have HOMER installed, the script can
# run annotatePeaks.pl on CHIP_UNION_BED and create this TSV.
CHIP_ANNOT_TSV = WORKDIR / "results/chromatin_example/union_annotation.tsv"

# Genome assembly used for annotation (e.g. "hg38", "mm10")
GENOME_ASSEMBLY = "hg38"

# If you want the script to auto-run HOMER annotatePeaks when CHIP_ANNOT_TSV
# is missing, set this to either:
#   - the directory containing annotatePeaks.pl
#   - or None to try auto-detect via environment modules (BioHPC style)
HOMER_BIN_DIR = None   # e.g. Path("/path/to/homer/bin") or None

# Column names inside the annotation table (CHIP_ANNOT_TSV).
# These defaults are compatible with HOMER annotatePeaks output.
ANNOT_CHR_COL      = "Chr"
ANNOT_START_COL    = "Start"
ANNOT_END_COL      = "End"
ANNOT_PEAKID_COL   = None       # if None, peak_id will be built as chr_start_end
ANNOT_SYMBOL_COL   = "Nearest PromoterID"
ANNOT_ENSEMBL_COL  = "Nearest Ensembl"
ANNOT_TYPE_COL     = "Annotation"  # used to identify "Promoter" peaks

# ----------------- RNA DESeq2 inputs ---------------------

# RNA DESeq2 result table (one file). If set to None, the script will
# choose the largest *.csv/tsv/txt file in RNADIR.
RNADIR = WORKDIR / "results/rna_example"
RNA_DE_FILE = None  # e.g. Path("/path/to/deseq2_rna_results.csv") or None

# Column names in RNA DE table
# It is safer to set them explicitly; if left as None, heuristic detection will be used.
RNA_ID_COL       = None   # e.g. "gene_id" / "ensembl" / "Geneid"
RNA_SYMBOL_COL   = None   # optional symbol column, e.g. "symbol" / "SYMBOL"
RNA_LOG2FC_COL   = None   # e.g. "log2FoldChange"
RNA_FDR_COL      = None   # e.g. "padj"

# ------------- Output directory and naming ----------------

OUTDIR  = WORKDIR / "results/integration_rna_chromatin"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Base filenames for outputs (no project-specific suffix by default)
SCATTER_PDF_NAME = "RNA_vs_chromatin_quadrants.pdf"
Q1_GENES_CSV     = "Quadrant_I_both_up_genes.csv"
Q3_GENES_CSV     = "Quadrant_III_both_down_genes.csv"
COUNTS_CSV       = "Quadrant_counts.csv"
JOINED_TABLE_CSV = "RNA_ChIP_promoter_joined.csv"

# ------ Statistics and visualization parameters ----------

FDR_THRESH   = 0.05
RNA_LFC_THR  = 1.0
CHIP_LFC_THR = 1.0

X_LIM   = (-6, 6)
Y_LIM   = (-6, 6)
X_LABEL = "Promoter chromatin log2FC (group2 / group1)"
Y_LABEL = "RNA log2FC (group2 / group1)"
TITLE   = "RNA vs chromatin (promoter level)"

# Optional gene marker lists (highlighted on scatter plot).
# Leave empty lists [] if you do not want to highlight anything.
MARKERS_SET_A = ["KRT5", "KRT14", "TP63"]   # e.g. basal markers
MARKERS_SET_B = ["KRT4", "KRT13", "KLF4"]   # e.g. luminal markers

# =========================================================
# END OF CONFIG
# =========================================================


def say(*args):
    print("[06_rna_chromatin]", *args)


# ---------------------------------------------------------
# 1) HOMER annotation: annotatePeaks on union BED (optional)
# ---------------------------------------------------------

def auto_detect_homer_bin():
    """
    Best-effort attempt to detect HOMER bin dir via environment modules.
    This is BioHPC-friendly but optional; if it fails, return None.
    """
    try:
        init = (
            "source /etc/profile.d/modules.sh 2>/dev/null || "
            "source /usr/share/Modules/init/bash 2>/dev/null || true"
        )
        cmd = (
            f"{init}; "
            "module load load_rhel7_apps || true; "
            "module load homer/4.9 || true; "
            "dirname $(which annotatePeaks.pl) 2>/dev/null || true"
        )
        res = subprocess.run(
            ["bash", "-lc", cmd], capture_output=True, text=True, check=False
        )
        guess = (res.stdout or "").strip()
        if guess and os.path.exists(guess):
            return Path(guess)
    except Exception:
        pass
    return None


def ensure_peak_annotation():
    """
    Ensure CHIP_ANNOT_TSV exists.
    If not, and HOMER_BIN_DIR is available (or auto-detected),
    run annotatePeaks.pl on CHIP_UNION_BED and save as TSV.
    """
    if CHIP_ANNOT_TSV.exists():
        say("Annotation TSV already exists:", CHIP_ANNOT_TSV)
        return

    if not CHIP_UNION_BED.exists():
        raise FileNotFoundError(f"Union BED not found: {CHIP_UNION_BED}")

    homer_bin = HOMER_BIN_DIR
    if homer_bin is None:
        homer_bin = auto_detect_homer_bin()

    if homer_bin is None:
        raise RuntimeError(
            "CHIP_ANNOT_TSV does not exist and HOMER_BIN_DIR is not set / not detectable.\n"
            "Please either provide an existing annotation TSV or set HOMER_BIN_DIR."
        )

    annotate = Path(homer_bin) / "annotatePeaks.pl"
    if not annotate.exists():
        raise FileNotFoundError(f"annotatePeaks.pl not found in {homer_bin}")

    out_txt = CHIP_ANNOT_TSV.with_suffix(".txt")
    cmd = (
        f"export PATH='{homer_bin}':$PATH; "
        f"{annotate} {CHIP_UNION_BED} {GENOME_ASSEMBLY} "
        f"-size given -cpu 8 > {out_txt}"
    )
    say("Running HOMER annotatePeaks.pl ...")
    say(cmd)
    subprocess.run(["bash", "-lc", cmd], check=False)

    say("Parsing HOMER annotation to TSV ...")
    df = pd.read_csv(out_txt, sep="\t", comment="#")
    # Add peak_id if missing, using chr_start_end
    if ANNOT_PEAKID_COL and ANNOT_PEAKID_COL in df.columns:
        df["peak_id"] = df[ANNOT_PEAKID_COL].astype(str)
    else:
        required = {ANNOT_CHR_COL, ANNOT_START_COL, ANNOT_END_COL}
        if not required.issubset(df.columns):
            raise ValueError(
                f"Annotation table missing required columns {required} "
                f"and no explicit peak_id column is specified."
            )
        df["peak_id"] = (
            df[ANNOT_CHR_COL].astype(str)
            + "_"
            + df[ANNOT_START_COL].astype(str)
            + "_"
            + df[ANNOT_END_COL].astype(str)
        )

    df.to_csv(CHIP_ANNOT_TSV, sep="\t", index=False)
    say("WROTE annotation TSV:", CHIP_ANNOT_TSV)


# ---------------------------------------------------------
# 2) Load RNA DE results
# ---------------------------------------------------------

def load_rna_de():
    """
    Load RNA DE result table and standardize columns:
    - rna_log2fc
    - rna_fdr
    - ensembl (if available)
    - rna_symbol (if available)
    """
    # Determine RNA file
    if RNA_DE_FILE is not None:
        rna_file = Path(RNA_DE_FILE)
        if not rna_file.exists():
            raise FileNotFoundError(f"RNA DE file not found: {rna_file}")
    else:
        candidates = (
            list(RNADIR.glob("*.csv")) +
            list(RNADIR.glob("*.tsv")) +
            list(RNADIR.glob("*.txt"))
        )
        if not candidates:
            raise FileNotFoundError(f"No RNA DE files found in {RNADIR}")
        rna_file = sorted(candidates, key=lambda p: p.stat().st_size, reverse=True)[0]
    sep = "\t" if rna_file.suffix.lower() in {".tsv", ".txt"} else ","
    say("RNA DE file:", rna_file)
    rna = pd.read_csv(rna_file, sep=sep)

    # Column detection (if not specified in CONFIG)
    log2fc_col = RNA_LOG2FC_COL
    fdr_col = RNA_FDR_COL
    id_col = RNA_ID_COL
    sym_col = RNA_SYMBOL_COL

    if log2fc_col is None:
        for cand in ["log2FoldChange", "log2FC", "log2fc"]:
            if cand in rna.columns:
                log2fc_col = cand
                break
    if fdr_col is None:
        for cand in ["padj", "FDR", "qvalue", "q_value"]:
            if cand in rna.columns:
                fdr_col = cand
                break
    if id_col is None:
        for cand in ["ensembl", "ENSG", "gene_id", "Geneid", "geneId", "gene_id_version"]:
            if cand in rna.columns:
                id_col = cand
                break
    if sym_col is None:
        for cand in ["symbol", "SYMBOL", "Gene", "gene", "GENE"]:
            if cand in rna.columns:
                sym_col = cand
                break

    if log2fc_col is None or fdr_col is None:
        raise ValueError(
            "Could not detect RNA log2FC / FDR columns. "
            "Please set RNA_LOG2FC_COL and RNA_FDR_COL in the CONFIG."
        )

    if id_col is None and sym_col is None:
        raise ValueError(
            "Could not detect any ID column (Ensembl or symbol) in RNA table. "
            "Please set RNA_ID_COL and/or RNA_SYMBOL_COL."
        )

    # Standardize column names
    rna = rna.rename(columns={log2fc_col: "rna_log2fc", fdr_col: "rna_fdr"})
    cols = ["rna_log2fc", "rna_fdr"]
    if id_col is not None:
        rna = rna.rename(columns={id_col: "ensembl"})
        rna["ensembl"] = rna["ensembl"].astype(str).str.strip()
        cols.append("ensembl")
    if sym_col is not None:
        rna = rna.rename(columns={sym_col: "rna_symbol"})
        rna["rna_symbol"] = rna["rna_symbol"].astype(str).str.strip()
        cols.append("rna_symbol")

    rna = rna[cols].dropna(subset=[c for c in cols if c in ["ensembl", "rna_symbol"]])
    say("Loaded RNA DE rows:", len(rna))
    return rna


# ---------------------------------------------------------
# 3) Load chromatin DE + annotation, aggregate to gene
# ---------------------------------------------------------

def load_chip_de():
    """
    Load chromatin DESeq2 results (peak-level) and merge with annotation.
    Returns a DataFrame with columns such as:
      peak_id, chip_log2fc, chip_fdr, symbol_label, ensembl, ann_type
    """
    # Choose DE file
    chip_file = None
    if CHIP_DE_UNSHRUNKEN.exists():
        chip_file = CHIP_DE_UNSHRUNKEN
    elif CHIP_DE_SHRUNKEN.exists():
        chip_file = CHIP_DE_SHRUNKEN
    else:
        # fallback: pick any CSV in directory
        candidates = list(CHIP_DE_DIR.glob("*.csv")) + list(CHIP_DE_DIR.glob("*.tsv"))
        if candidates:
            chip_file = sorted(candidates, key=lambda p: p.stat().st_size, reverse=True)[0]
        else:
            raise FileNotFoundError(f"No chromatin DE files found in {CHIP_DE_DIR}")

    sep = "\t" if chip_file.suffix.lower() in {".tsv", ".txt"} else ","
    say("Chromatin DE file:", chip_file)
    chip = pd.read_csv(chip_file, sep=sep)

    # log2FC / FDR columns
    if "log2FoldChange_raw" in chip.columns:
        lfc_col = "log2FoldChange_raw"
    else:
        lfc_col = None
        for cand in ["log2FoldChange", "log2FC", "log2fc"]:
            if cand in chip.columns:
                lfc_col = cand
                break
    fdr_col = None
    for cand in ["padj_raw", "padj", "FDR", "qvalue", "q_value"]:
        if cand in chip.columns:
            fdr_col = cand
            break

    if "peak_id" not in chip.columns or lfc_col is None or fdr_col is None:
        raise ValueError(
            "Chromatin DE table must contain 'peak_id' and log2FC/FDR columns. "
            "Detected: peak_id in %s, lfc_col=%s, fdr_col=%s"
            % (chip.columns, lfc_col, fdr_col)
        )

    chip = chip.rename(columns={lfc_col: "chip_log2fc", fdr_col: "chip_fdr"})

    # Annotation
    if not CHIP_ANNOT_TSV.exists():
        raise FileNotFoundError(
            f"Annotation TSV not found: {CHIP_ANNOT_TSV}. "
            "Run ensure_peak_annotation() first."
        )

    annot = pd.read_csv(CHIP_ANNOT_TSV, sep="\t")
    # Build peak_id if needed
    if "peak_id" not in annot.columns:
        required = {ANNOT_CHR_COL, ANNOT_START_COL, ANNOT_END_COL}
        if not required.issubset(annot.columns):
            raise ValueError(
                f"Annotation table missing columns {required} and no 'peak_id' column."
            )
        annot["peak_id"] = (
            annot[ANNOT_CHR_COL].astype(str)
            + "_"
            + annot[ANNOT_START_COL].astype(str)
            + "_"
            + annot[ANNOT_END_COL].astype(str)
        )

    # Symbol / Ensembl / annotation type
    sym_series = (
        annot[ANNOT_SYMBOL_COL].astype(str)
        if ANNOT_SYMBOL_COL in annot.columns
        else pd.Series([""] * len(annot))
    )
    ens_series = (
        annot[ANNOT_ENSEMBL_COL].astype(str)
        if ANNOT_ENSEMBL_COL in annot.columns
        else pd.Series([""] * len(annot))
    )
    ann_type_series = (
        annot[ANNOT_TYPE_COL].astype(str)
        if ANNOT_TYPE_COL in annot.columns
        else pd.Series([""] * len(annot))
    )

    a = annot.copy()
    a["symbol_label"] = (
        sym_series.fillna("")
        .str.split("[,;/\\s()]", regex=True)
        .str[0]
        .str.replace("^NA$", "", regex=True)
        .str.strip()
    )
    a["ensembl"] = (
        ens_series.fillna("")
        .str.split("[,;/\\s()]", regex=True)
        .str[0]
        .str.replace("^NA$", "", regex=True)
        .str.strip()
    )
    a["ann_type"] = ann_type_series

    merge_cols = ["peak_id", "symbol_label", "ensembl", "ann_type"]
    annot_small = a[merge_cols].copy()

    chip = chip.merge(annot_small, on="peak_id", how="left")
    say("Chromatin DE rows:", len(chip))
    return chip


def aggregate_chip_to_gene(chip_df):
    """
    Aggregate peak-level chromatin DE to gene level:
    - best peak per gene (min FDR) -> chip_log2fc_best, chip_fdr_best
    - promoter-only peak if available (max |LFC| promoter) -> chip_log2fc_prom, chip_fdr_prom
    - also provide symbol_label lookup per Ensembl
    """
    # Best peak (min FDR)
    chip_best = (
        chip_df.dropna(subset=["ensembl"])
        .sort_values(["ensembl", "chip_fdr"], ascending=[True, True])
        .groupby("ensembl", as_index=False)
        .first()[["ensembl", "chip_log2fc", "chip_fdr"]]
        .rename(columns={"chip_log2fc": "chip_log2fc_best", "chip_fdr": "chip_fdr_best"})
    )

    # Promoter-only aggregation if annotation available
    if "ann_type" in chip_df.columns:
        prom = chip_df[
            chip_df["ann_type"].astype(str).str.contains("Promoter", case=False, na=False)
        ].copy()
        if len(prom):
            prom["_abs"] = prom["chip_log2fc"].abs()
            chip_prom = (
                prom.dropna(subset=["ensembl"])
                .sort_values(["ensembl", "_abs"], ascending=[True, False])
                .groupby("ensembl", as_index=False)
                .first()[["ensembl", "chip_log2fc", "chip_fdr"]]
                .rename(columns={"chip_log2fc": "chip_log2fc_prom", "chip_fdr": "chip_fdr_prom"})
            )
        else:
            chip_prom = pd.DataFrame(columns=["ensembl", "chip_log2fc_prom", "chip_fdr_prom"])
    else:
        chip_prom = pd.DataFrame(columns=["ensembl", "chip_log2fc_prom", "chip_fdr_prom"])

    # Symbol lookup
    symbol_lookup = (
        chip_df[["ensembl", "symbol_label"]]
        .dropna()
        .drop_duplicates(subset=["ensembl"])
    )

    chip_gene = chip_best.merge(chip_prom, on="ensembl", how="left")
    chip_gene = chip_gene.merge(symbol_lookup, on="ensembl", how="left")
    say("Gene-level chromatin rows (Ensembl):", len(chip_gene))
    return chip_gene


# ---------------------------------------------------------
# 4) Join RNA & chromatin, classify quadrants, plot
# ---------------------------------------------------------

def join_and_plot(rna_df, chip_gene_df):
    """
    Join RNA and chromatin at gene-level, classify quadrants, and plot.
    """
    # Use promoter-only log2FC if available; else fall back to best-peak.
    if {"chip_log2fc_prom", "chip_fdr_prom"}.issubset(chip_gene_df.columns):
        chip_prom = chip_gene_df[
            ["ensembl", "symbol_label", "chip_log2fc_prom", "chip_fdr_prom"]
        ].dropna(subset=["ensembl"])
        chip_prom = chip_prom.rename(
            columns={"chip_log2fc_prom": "chip_log2fc", "chip_fdr_prom": "chip_fdr"}
        )
        say("Using promoter-only aggregation for chromatin.")
    else:
        chip_prom = chip_gene_df[
            ["ensembl", "symbol_label", "chip_log2fc_best", "chip_fdr_best"]
        ].dropna(subset=["ensembl"])
        chip_prom = chip_prom.rename(
            columns={"chip_log2fc_best": "chip_log2fc", "chip_fdr_best": "chip_fdr"}
        )
        say("Promoter-only aggregation not available; using best peak per gene.")

    # Join on Ensembl
    if "ensembl" not in rna_df.columns:
        raise ValueError("RNA table has no 'ensembl' column; cannot perform Ensembl-based join.")

    joined = rna_df.merge(chip_prom, on="ensembl", how="inner").dropna()
    say("Joined RNA–chromatin rows:", len(joined))
    if joined.empty:
        say("[WARN] No overlapping genes after join. Skipping scatter plot and outputs.")
        return

    # Spearman correlation
    if HAS_SCIPY:
        mask = np.isfinite(joined["rna_log2fc"]) & np.isfinite(joined["chip_log2fc"])
        rho, pval = spearmanr(joined.loc[mask, "rna_log2fc"], joined.loc[mask, "chip_log2fc"])
        say(f"Spearman rho={rho:.3f}, p={pval:.2e}")
    else:
        rho, pval = np.nan, np.nan
        say("scipy not available; Spearman correlation skipped.")

    # Quadrant classification
    cond_up_up = (
        (joined["rna_fdr"] <= FDR_THRESH) &
        (joined["chip_fdr"] <= FDR_THRESH) &
        (joined["rna_log2fc"] >= RNA_LFC_THR) &
        (joined["chip_log2fc"] >= CHIP_LFC_THR)
    )
    cond_dn_dn = (
        (joined["rna_fdr"] <= FDR_THRESH) &
        (joined["chip_fdr"] <= FDR_THRESH) &
        (joined["rna_log2fc"] <= -RNA_LFC_THR) &
        (joined["chip_log2fc"] <= -CHIP_LFC_THR)
    )
    cond_oppo = (
        (joined["rna_fdr"] <= FDR_THRESH) &
        (joined["chip_fdr"] <= FDR_THRESH) &
        (
            ((joined["rna_log2fc"] >= RNA_LFC_THR) & (joined["chip_log2fc"] <= -CHIP_LFC_THR)) |
            ((joined["rna_log2fc"] <= -RNA_LFC_THR) & (joined["chip_log2fc"] >= CHIP_LFC_THR))
        )
    )
    state = np.where(
        cond_up_up, "Both_Up",
        np.where(cond_dn_dn, "Both_Down",
                 np.where(cond_oppo, "Opposite", "NS"))
    )
    joined["quadrant"] = state

    q_counts = {
        "Both_Up": int(cond_up_up.sum()),
        "Both_Down": int(cond_dn_dn.sum()),
        "Opposite": int(cond_oppo.sum()),
        "NS": int((joined["quadrant"] == "NS").sum()),
        "Total": int(len(joined)),
    }
    say("Quadrant counts:", q_counts)

    # Scatter plot
    plt.figure(figsize=(6.5, 6))
    # background
    plt.scatter(
        joined["chip_log2fc"], joined["rna_log2fc"],
        c="#bfbfbf", s=8, alpha=0.25, linewidths=0, zorder=1
    )
    # Quadrant I (both up)
    plt.scatter(
        joined.loc[cond_up_up, "chip_log2fc"],
        joined.loc[cond_up_up, "rna_log2fc"],
        c="#d62728", s=12, alpha=0.95, linewidths=0, zorder=2,
        label="Quadrant I (both up)"
    )
    # Quadrant III (both down)
    plt.scatter(
        joined.loc[cond_dn_dn, "chip_log2fc"],
        joined.loc[cond_dn_dn, "rna_log2fc"],
        c="#1f77b4", s=12, alpha=0.95, linewidths=0, zorder=2,
        label="Quadrant III (both down)"
    )

    # Threshold lines
    plt.axvline(CHIP_LFC_THR, color="gray", ls=(0, (6, 4)), lw=1.2)
    plt.axvline(-CHIP_LFC_THR, color="gray", ls=(0, (6, 4)), lw=1.2)
    plt.axhline(RNA_LFC_THR, color="gray", ls=(0, (6, 4)), lw=1.2)
    plt.axhline(-RNA_LFC_THR, color="gray", ls=(0, (6, 4)), lw=1.2)

    plt.xlim(*X_LIM)
    plt.ylim(*Y_LIM)
    plt.xlabel(X_LABEL)
    plt.ylabel(Y_LABEL)
    if np.isfinite(rho):
        plt.title(f"{TITLE}\nSpearman rho={rho:.3f}, p={pval:.2e}")
    else:
        plt.title(TITLE)

    # Marker highlighting (optional)
    label_upper = joined["symbol_label"].astype(str).str.upper() if "symbol_label" in joined.columns else pd.Series([""] * len(joined))

    # Build symbol → ensembl map from chromatin gene table
    sym2ens_map = {}
    if "symbol_label" in chip_gene_df.columns:
        tmp = (
            chip_gene_df[["ensembl", "symbol_label"]]
            .dropna()
            .assign(SYM=lambda d: d["symbol_label"].astype(str).str.upper().str.strip())
            .drop_duplicates(subset=["SYM"])
        )
        sym2ens_map = tmp.set_index("SYM")["ensembl"].to_dict()

    def pick_one(sym, cond_mask):
        """
        Choose one index for a given gene symbol and quadrant mask:
        1) try Ensembl-based mapping (more robust),
        2) fall back to symbol_label exact / contains match.
        """
        sym_up = sym.upper()
        # 1) via Ensembl mapping
        ensg = sym2ens_map.get(sym_up)
        if ensg is not None and "ensembl" in joined.columns:
            idxs = np.where(joined["ensembl"].astype(str) == ensg)[0]
            if len(idxs):
                cand = idxs[cond_mask.iloc[idxs]]
                if len(cand) == 0:
                    cand = idxs
                j = cand[np.nanargmin(joined["rna_fdr"].iloc[cand].values)]
                return int(j)
        # 2) fallback: match symbol_label
        idxs = np.where(label_upper == sym_up)[0]
        if len(idxs) == 0:
            idxs = np.where(label_upper.str.contains(rf"\b{sym_up}\b", regex=True))[0]
        if len(idxs) == 0:
            return None
        cand = idxs[cond_mask.iloc[idxs]]
        if len(cand) == 0:
            cand = idxs
        j = cand[np.nanargmin(joined["rna_fdr"].iloc[cand].values)]
        return int(j)

    labeled = []

    # Markers in Quadrant I (both up), e.g. luminal markers
    for sym in MARKERS_SET_B:
        j = pick_one(sym, cond_up_up)
        if j is None:
            continue
        x = float(joined.at[j, "chip_log2fc"])
        y = float(joined.at[j, "rna_log2fc"])
        txt = plt.annotate(
            sym,
            (x, y),
            xytext=(x + 0.15, y + 0.15),
            textcoords="data",
            color="black",
            zorder=5,
            arrowprops=dict(
                arrowstyle="-", lw=0.9, color="black", shrinkA=0, shrinkB=0
            ),
            fontsize=10,
        )
        txt.set_path_effects([pe.withStroke(linewidth=2, foreground="white")])
        labeled.append(sym)

    # Markers in Quadrant III (both down), e.g. basal markers
    for sym in MARKERS_SET_A:
        j = pick_one(sym, cond_dn_dn)
        if j is None:
            continue
        x = float(joined.at[j, "chip_log2fc"])
        y = float(joined.at[j, "rna_log2fc"])
        txt = plt.annotate(
            sym,
            (x, y),
            xytext=(x - 0.15, y - 0.15),
            textcoords="data",
            color="black",
            zorder=5,
            arrowprops=dict(
                arrowstyle="-", lw=0.9, color="black", shrinkA=0, shrinkB=0
            ),
            fontsize=10,
        )
        txt.set_path_effects([pe.withStroke(linewidth=2, foreground="white")])
        labeled.append(sym)

    say("Labeled markers:", labeled)

    plt.tight_layout()
    scatter_pdf = OUTDIR / SCATTER_PDF_NAME
    plt.savefig(scatter_pdf)
    plt.close()
    say("WROTE scatter plot:", scatter_pdf)

    # Export quadrant gene tables
    q1 = joined.loc[
        cond_up_up,
        ["ensembl", "symbol_label", "rna_log2fc", "rna_fdr", "chip_log2fc", "chip_fdr"],
    ].copy()
    q3 = joined.loc[
        cond_dn_dn,
        ["ensembl", "symbol_label", "rna_log2fc", "rna_fdr", "chip_log2fc", "chip_fdr"],
    ].copy()

    q1_out = OUTDIR / Q1_GENES_CSV
    q3_out = OUTDIR / Q3_GENES_CSV
    q1.to_csv(q1_out, index=False)
    q3.to_csv(q3_out, index=False)
    say("WROTE QI genes:", q1_out)
    say("WROTE QIII genes:", q3_out)

    # Export quadrant counts
    counts_out = OUTDIR / COUNTS_CSV
    pd.DataFrame(
        [
            {"quadrant": "Quadrant I (Both_Up)", "n": int(cond_up_up.sum())},
            {"quadrant": "Quadrant III (Both_Down)", "n": int(cond_dn_dn.sum())},
            {"quadrant": "Opposite", "n": int(cond_oppo.sum())},
            {"quadrant": "NS", "n": int((joined["quadrant"] == "NS").sum())},
            {"quadrant": "Total", "n": int(len(joined))},
        ]
    ).to_csv(counts_out, index=False)
    say("WROTE quadrant counts:", counts_out)

    # Export full joined table
    joined_out = OUTDIR / JOINED_TABLE_CSV
    joined.to_csv(joined_out, index=False)
    say("WROTE joined table:", joined_out)


# ---------------------------------------------------------
# main
# ---------------------------------------------------------

if __name__ == "__main__":
    # 1) Ensure annotation
    ensure_peak_annotation()

    # 2) Load RNA & chromatin
    rna_df = load_rna_de()
    chip_df = load_chip_de()
    chip_gene_df = aggregate_chip_to_gene(chip_df)

    # 3) Join + plot
    join_and_plot(rna_df, chip_gene_df)

    say("DONE.")
