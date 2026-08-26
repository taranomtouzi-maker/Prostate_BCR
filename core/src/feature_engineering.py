"""
Feature engineering for Prostate BCR prediction.

Extracts clinically meaningful interaction and pathway features.
This module is the single source of truth for all engineered features.

All feature engineering must be fitted on training data only.
When used in external validation, pathway scores are computed on the
common gene space to ensure transferability.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────
# Default gene sets (import from pathways.py for consistency)
# ──────────────────────────────────────────────────────────────────────

def _default_gene_sets() -> Dict[str, List[str]]:
    """Import pathway gene sets lazily to avoid circular imports."""
    from src.pathways import PSA_PATHWAY, AR_SIGNALING, PROLIFERATION
    return {
        "PSA": PSA_PATHWAY,
        "AR": AR_SIGNALING,
        "PROLIF": PROLIFERATION,
    }


# ──────────────────────────────────────────────────────────────────────
# Clinical column constants
# ──────────────────────────────────────────────────────────────────────
GLEASON_PRIMARY_COL = "Gleason pattern primary"
GLEASON_SECONDARY_COL = "Gleason pattern secondary"
MARGIN_COL = "Surgical Margin Resection Status_R1"
LYMPH_NODE_COL = "Primary Lymph Node Presentation Assessment Ind-3_YES"


# ──────────────────────────────────────────────────────────────────────
# Engineered feature creation
# ──────────────────────────────────────────────────────────────────────

def create_engineered_features(
    X: pd.DataFrame,
    *,
    strict_mode: bool = False,
    gene_sets: Dict[str, List[str]] | None = None,
    min_genes_for_pathway: int = 3,
) -> Tuple[pd.DataFrame, List[str]]:
    """Create clinically meaningful engineered features.

    Creates:
      1. Gleason_Total (primary + secondary)
      2. High_Risk_Gleason (any pattern >= 4)
      3. Margin_x_LymphNode interaction
      4. T_Stage_Risk (sum of T3/T4 dummy variables)
      5. PSA_Pathway_Score (mean of PSA genes)
      6. AR_Signaling_Score (mean of AR pathway genes)
      7. Proliferation_Score (mean of proliferation genes)

    Parameters
    ----------
    X : Input DataFrame.
    strict_mode : If True, only create pathway scores if ALL genes present.
    gene_sets : Optional dict overriding default gene sets.
    min_genes_for_pathway : Minimum genes required for a pathway score.

    Returns
    -------
    Tuple of (DataFrame with new features, list of new feature names).
    """
    X = X.copy()
    created_features: List[str] = []
    gene_sets = gene_sets or _default_gene_sets()

    # ── 1. Gleason-based features ──
    if GLEASON_PRIMARY_COL in X.columns and GLEASON_SECONDARY_COL in X.columns:
        X["Gleason_Total"] = X[GLEASON_PRIMARY_COL] + X[GLEASON_SECONDARY_COL]
        X["High_Risk_Gleason"] = (
            (X[GLEASON_PRIMARY_COL] >= 4) | (X[GLEASON_SECONDARY_COL] >= 4)
        ).astype(int)
        created_features.extend(["Gleason_Total", "High_Risk_Gleason"])
        logger.info("Engineered: Gleason_Total, High_Risk_Gleason")

    # ── 2. Margin × Lymph Node interaction ──
    if MARGIN_COL in X.columns and LYMPH_NODE_COL in X.columns:
        X["Margin_x_LymphNode"] = (
            X[MARGIN_COL].astype(float) * X[LYMPH_NODE_COL].astype(float)
        )
        created_features.append("Margin_x_LymphNode")
        logger.info("Engineered: Margin_x_LymphNode")

    # ── 3. T-Stage risk score ──
    t_stage_cols = [c for c in X.columns if "Tumor Stage Code_T3" in c or "Tumor Stage Code_T4" in c]
    if len(t_stage_cols) >= 2:
        X["T_Stage_Risk"] = X[t_stage_cols].sum(axis=1)
        created_features.append("T_Stage_Risk")
        logger.info("Engineered: T_Stage_Risk")

    # ── 4–6. Pathway scores ──
    pathway_names = {"PSA": "PSA_Pathway_Score", "AR": "AR_Signaling_Score", "PROLIF": "Proliferation_Score"}
    for key, score_name in pathway_names.items():
        genes = gene_sets.get(key, [])
        available = [g for g in genes if g in X.columns]

        if strict_mode:
            if set(genes).issubset(set(X.columns)):
                X[score_name] = X[genes].mean(axis=1)
                created_features.append(score_name)
                logger.info("Engineered: %s (strict, %d genes)", score_name, len(genes))
        else:
            if len(available) >= min_genes_for_pathway:
                X[score_name] = X[available].mean(axis=1)
                created_features.append(score_name)
                logger.info("Engineered: %s (from %d genes)", score_name, len(available))

    return X, created_features


def create_pathway_score_features(
    X: pd.DataFrame,
    pathway_scores: pd.DataFrame,
) -> Tuple[pd.DataFrame, List[str]]:
    """Append pathway score columns to a feature DataFrame.

    Parameters
    ----------
    X : Base feature DataFrame.
    pathway_scores : DataFrame with pathway score columns (same index as X).

    Returns
    -------
    Tuple of (merged DataFrame, list of added pathway column names).
    """
    added_cols = [c for c in pathway_scores.columns if c not in X.columns]
    if not added_cols:
        return X, []

    X = X.copy()
    for col in added_cols:
        X[col] = pathway_scores[col]

    logger.info("Added %d pathway scores: %s", len(added_cols), added_cols)
    return X, added_cols


# ──────────────────────────────────────────────────────────────────────
# Validation helpers
# ──────────────────────────────────────────────────────────────────────

def validate_engineered_features(
    X_train: pd.DataFrame,
    X_external: pd.DataFrame,
    engineered_cols: List[str],
) -> Dict[str, any]:
    """Check which engineered features are available in external data.

    Returns a dict with available/missing/transferable counts.
    """
    available_in_ext = [c for c in engineered_cols if c in X_external.columns]
    missing_in_ext = [c for c in engineered_cols if c not in X_external.columns]

    # Separate clinical vs gene-based engineered features
    clinical_engineered = ["Gleason_Total", "High_Risk_Gleason", "Margin_x_LymphNode", "T_Stage_Risk"]
    pathway_engineered = ["PSA_Pathway_Score", "AR_Signaling_Score", "Proliferation_Score"]

    result = {
        "total": len(engineered_cols),
        "available_in_external": len(available_in_ext),
        "missing_in_external": missing_in_ext,
        "clinical_available": [c for c in clinical_engineered if c in available_in_ext],
        "pathway_available": [c for c in pathway_engineered if c in available_in_ext],
    }

    if missing_in_ext:
        logger.warning(
            "External validation: %d engineered features missing: %s",
            len(missing_in_ext), missing_in_ext,
        )

    return result
