"""
Data loaders for external validation cohorts.

Provides robust loading functions for GEO datasets (GSE70769, GSE54460)
with:
  - Correct BCR label mapping (y/n/n/a)
  - Missing clinical column handling with explicit warnings
  - Cross-platform normalization (quantile, patient-zscore, frozen ComBat)
  - Feature intersection with training data

All functions are designed for use in the consolidated evaluation notebook.
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.preprocessing import quantile_transform

import config
from src.io import logger

# ──────────────────────────────────────────────────────────────────────
# GSE70769 Loader
# ──────────────────────────────────────────────────────────────────────

# Expected clinical columns in GSE70769 (after parsing)
# These may be missing → we handle gracefully with warnings
EXPECTED_CLINICAL_COLS = [
    "Margin",             # Surgical margin status
    "Lymph_Node",         # Lymph node involvement
    "Gleason",            # Gleason score
    "PSA",                # Pre-operative PSA
    "Age",                # Patient age
    "T_Stage",            # Tumor stage
]


def load_gse70769(
    base_dir: Path | str | None = None,
    *,
    gene_expression_path: Path | str | None = None,
    clinical_path: Path | str | None = None,
    bcr_col: str = "BCR",
    normalize: str = "quantile",
    common_genes: Sequence[str] | None = None,
) -> Tuple[pd.DataFrame, pd.Series, Dict]:
    """Load and preprocess GSE70769 external validation cohort.

    Handles:
      - BCR label mapping ('y'→1, 'n'→0, 'n/a'→NaN→drop)
      - Missing clinical columns with explicit warnings
      - Cross-platform normalization (quantile or rank)
      - Feature intersection with training gene space

    Parameters
    ----------
    base_dir : Root directory for GSE70769 data files.
    gene_expression_path : Path to gene expression CSV/TSV.
    clinical_path : Path to clinical metadata CSV/TSV.
    bcr_col : Column name containing BCR status.
    normalize : Normalization method ('quantile', 'rank', 'zscore', 'none').
    common_genes : If provided, restrict to intersection with these genes.

    Returns
    -------
    Tuple of (X_external, y_external, metadata_dict).
    """
    base_dir = Path(base_dir or config.PROCESSED_DIR)

    # ── Load gene expression ──
    if gene_expression_path is not None:
        expr_path = Path(gene_expression_path)
    else:
        # Try common locations
        candidates = [
            base_dir / "X_GSE70769.csv",
            base_dir / "GSE70769_expression.csv",
            config.DATA_DIR / "external" / "X_GSE70769.csv",
        ]
        expr_path = None
        for c in candidates:
            if c.exists():
                expr_path = c
                break
        if expr_path is None:
            raise FileNotFoundError(
                f"GSE70769 expression file not found. Tried: {[str(c) for c in candidates]}"
            )

    logger.info("Loading GSE70769 expression from %s", expr_path)

    if expr_path.suffix == ".tsv":
        X_expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    else:
        X_expr = pd.read_csv(expr_path, index_col=0)

    logger.info("GSE70769 expression: %d samples × %d genes", *X_expr.shape)

    # ── Load clinical / labels ──
    if clinical_path is not None:
        clin_path = Path(clinical_path)
    else:
        candidates = [
            base_dir / "y_GSE70769.csv",
            base_dir / "GSE70769_clinical.csv",
            config.DATA_DIR / "external" / "y_GSE70769.csv",
        ]
        clin_path = None
        for c in candidates:
            if c.exists():
                clin_path = c
                break
        if clin_path is None:
            # Try to use expression index as patient IDs and look for BCR column
            logger.warning("No clinical file found; attempting to extract labels from expression index")
            clin_path = None

    metadata: Dict = {"n_raw_samples": len(X_expr)}

    if clin_path is not None:
        logger.info("Loading GSE70769 clinical from %s", clin_path)
        if clin_path.suffix == ".tsv":
            y_df = pd.read_csv(clin_path, sep="\t", index_col=0)
        else:
            y_df = pd.read_csv(clin_path, index_col=0)

        # ── Map BCR labels ──
        if bcr_col in y_df.columns:
            raw_labels = y_df[bcr_col]
        else:
            # Try to find a BCR-like column
            bcr_candidates = [c for c in y_df.columns if "bcr" in c.lower() or "recurrence" in c.lower()]
            if bcr_candidates:
                bcr_col = bcr_candidates[0]
                raw_labels = y_df[bcr_col]
                logger.info("Using column '%s' as BCR indicator", bcr_col)
            else:
                raise ValueError(f"No BCR column found in {clin_path}. Available: {list(y_df.columns[:10])}")

        # Map y/n/Y/N/True/False/1/0 → binary
        label_map = {
            "y": 1, "Y": 1, "yes": 1, "YES": 1, "Yes": 1,
            "True": 1, "TRUE": 1, "true": 1,
            "1": 1, 1: 1,
            "n": 0, "N": 0, "no": 0, "NO": 0, "No": 0,
            "False": 0, "FALSE": 0, "false": 0,
            "0": 0, 0: 0,
        }

        y = raw_labels.map(label_map)

        # Handle n/a, NaN, etc.
        n_before = len(y)
        n_invalid = y.isna().sum()
        if n_invalid > 0:
            # Check if original values were n/a-like
            raw_str = raw_labels.astype(str).str.lower()
            n_na = raw_str.isin(["n/a", "na", "nan", "none", "", "unk", "unknown"]).sum()
            logger.warning(
                "GSE70769: %d samples with unmapable BCR labels (n/a or missing) — dropping",
                max(n_invalid, n_na),
            )
            valid_mask = y.notna()
            y = y[valid_mask]
            X_expr = X_expr.loc[y.index]

        y = y.astype(int)
        metadata["n_after_label_filter"] = len(y)
        metadata["bcr_positive"] = int(y.sum())
        metadata["bcr_negative"] = int((y == 0).sum())
        metadata["bcr_rate"] = float(y.mean())

        logger.info(
            "GSE70769 labels: %d BCR+ / %d BCR− (%.1f%% positive rate)",
            metadata["bcr_positive"], metadata["bcr_negative"],
            metadata["bcr_rate"] * 100,
        )

        # ── Check clinical columns ──
        if "Margin" not in y_df.columns:
            logger.warning("GSE70769: 'Margin' column not found — will be unavailable for clinical model")
            metadata["missing_clinical"] = ["Margin"]
        else:
            metadata["missing_clinical"] = []

        # Check other expected clinical columns
        for col in EXPECTED_CLINICAL_COLS:
            if col not in y_df.columns and col != "Margin":
                metadata.setdefault("missing_clinical", []).append(col)
                logger.warning("GSE70769: clinical column '%s' not found", col)

        if metadata.get("missing_clinical"):
            logger.warning(
                "GSE70769 missing clinical columns: %s — clinical-only model will not be available",
                metadata["missing_clinical"],
            )
    else:
        # No clinical file; assume labels are in expression index or need external provision
        raise ValueError(
            "Clinical/label file required for GSE70769. "
            "Provide clinical_path or place y_GSE70769.csv in data/external/."
        )

    # ── Feature intersection ──
    if common_genes is not None:
        common = [g for g in common_genes if g in X_expr.columns]
        logger.info("GSE70769: %d / %d common genes with training data", len(common), len(common_genes))
        X_expr = X_expr[common]
        metadata["n_common_genes"] = len(common)
    else:
        metadata["n_common_genes"] = X_expr.shape[1]

    # ── Cross-platform normalization ──
    if normalize == "quantile":
        X_expr = _quantile_normalize_cohort(X_expr)
        logger.info("Applied quantile normalization to GSE70769")
    elif normalize == "rank":
        X_expr = _rank_normalize_cohort(X_expr)
        logger.info("Applied rank normalization to GSE70769")
    elif normalize == "zscore":
        X_expr = (X_expr - X_expr.mean()) / (X_expr.std() + 1e-8)
        logger.info("Applied z-score normalization to GSE70769")
    elif normalize != "none":
        logger.warning("Unknown normalization method '%s'; using raw values", normalize)

    metadata["normalization"] = normalize

    return X_expr, y, metadata


# ──────────────────────────────────────────────────────────────────────
# GSE54460 Loader (triangulation cohort)
# ──────────────────────────────────────────────────────────────────────

def load_gse54460(
    base_dir: Path | str | None = None,
    *,
    gene_expression_path: Path | str | None = None,
    clinical_path: Path | str | None = None,
) -> Tuple[pd.DataFrame, pd.Series, Dict]:
    """Load GSE54460 (RNA-Seq FFPE) for triangulation validation.

    Returns
    -------
    Tuple of (X_54460, y_54460, metadata_dict).
    """
    base_dir = Path(base_dir or config.DATA_DIR / "interim")

    expr_candidates = [
        Path(gene_expression_path) if gene_expression_path else None,
        base_dir / "X_GSE54460.csv",
        config.DATA_DIR / "interim" / "X_GSE54460.csv",
    ]
    expr_path = next((p for p in expr_candidates if p is not None and p.exists()), None)
    if expr_path is None:
        raise FileNotFoundError("GSE54460 expression file not found")

    clin_candidates = [
        Path(clinical_path) if clinical_path else None,
        base_dir / "y_GSE54460.csv",
        config.DATA_DIR / "interim" / "y_GSE54460.csv",
    ]
    clin_path = next((p for p in clin_candidates if p is not None and p.exists()), None)
    if clin_path is None:
        raise FileNotFoundError("GSE54460 clinical/label file not found")

    X = pd.read_csv(expr_path, index_col=0)
    y = pd.read_csv(clin_path, index_col=0).iloc[:, 0]

    metadata = {
        "n_samples": len(y),
        "n_genes": X.shape[1],
        "bcr_positive": int(y.sum()),
        "bcr_rate": float(y.mean()),
    }

    logger.info("GSE54460: %d samples, %d genes, BCR rate=%.1f%%",
                metadata["n_samples"], metadata["n_genes"], metadata["bcr_rate"] * 100)

    return X, y, metadata


# ──────────────────────────────────────────────────────────────────────
# Cross-platform normalization helpers
# ──────────────────────────────────────────────────────────────────────

def _quantile_normalize_cohort(X: pd.DataFrame) -> pd.DataFrame:
    """Apply quantile normalization (per-gene to normal distribution)."""
    values = X.values.astype(np.float64)
    normalized = quantile_transform(
        values, output_distribution="normal",
        subsample=100000, random_state=config.RANDOM_STATE,
    )
    return pd.DataFrame(normalized, index=X.index, columns=X.columns)


def _rank_normalize_cohort(X: pd.DataFrame) -> pd.DataFrame:
    """Convert expression to percentile ranks (0–1) per gene."""
    ranked = X.rank(method="average")
    return (ranked - ranked.min()) / (ranked.max() - ranked.min() + 1e-8)


def patient_zscore(
    df: pd.DataFrame,
    genes: Sequence[str],
) -> pd.DataFrame:
    """Per-patient z-score across given genes.

    Row-wise standardization removes per-array / per-library scale effects
    between RNA-Seq and microarray, leaving only relative expression.
    """
    sub = df[list(genes)].astype(float)
    mu = sub.mean(axis=1)
    sd = sub.std(axis=1, ddof=0).replace(0, 1e-9)
    return sub.sub(mu, axis=0).div(sd, axis=0)


def frozen_combat(
    source: pd.DataFrame,
    target: pd.DataFrame,
    genes: Sequence[str],
) -> pd.DataFrame:
    """Frozen ComBat: estimate location/scale on source, apply to target.

    Standardizes target with its OWN batch statistics, then re-applies the
    source (training) parameters — the standard frozen ComBat transfer.
    No target labels are touched, so validation stays honest.
    """
    src = source[list(genes)].astype(float)
    tgt = target[list(genes)].astype(float)

    src_mu = src.mean(axis=0)
    src_sd = src.std(axis=0, ddof=1).replace(0, 1e-9)
    tgt_mu = tgt.mean(axis=0)
    tgt_sd = tgt.std(axis=0, ddof=1).replace(0, 1e-9)

    standardized = tgt.sub(tgt_mu, axis=1).div(tgt_sd, axis=1)
    return standardized.mul(src_sd, axis=1).add(src_mu, axis=1)


def common_gene_space(
    train_cols: Sequence[str],
    external_cols: Sequence[str],
) -> List[str]:
    """Return intersection of gene columns across platforms."""
    external_set = set(external_cols)
    return [c for c in train_cols if c in external_set]


def get_common_features(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    exclude_patterns: Optional[List[str]] = None,
) -> List[str]:
    """Get common features between two DataFrames, optionally excluding patterns."""
    common = set(df1.columns) & set(df2.columns)
    if exclude_patterns:
        for pattern in exclude_patterns:
            common = {c for c in common if pattern not in c}
    return sorted(common)
