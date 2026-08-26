"""
Improved BCR prediction pipeline (implements IMPROVEMENT_PLAN.md, priorities 1-2).

Key corrections over the current pipeline:
  1. Stability selection (bootstrap L1) instead of PSO wrapper selection.
  2. Honest nested CV: feature selection runs INSIDE each outer fold.
  3. Calibrated Elastic Net as primary model (standard for p >> n genomics).
  4. Clinical-only baseline as the reference any genomic model must beat.
  5. Cross-platform external validation: common-feature-space training,
     per-patient z-score normalization, frozen ComBat, prevalence-aware
     threshold recalibration.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.utils import resample

logger = logging.getLogger(__name__)

RANDOM_STATE = 42

# ---------------------------------------------------------------------------
# Clinical baseline features (all RP-available, platform-independent)
# ---------------------------------------------------------------------------
CLINICAL_BASELINE_FEATURES: List[str] = [
    "Gleason pattern primary",
    "Gleason pattern secondary",
    "Radical Prostatectomy Gleason Score for Prostate Cancer",
    "Surgical Margin Resection Status_R1",
    "Neoplasm American Joint Committee on Cancer Clinical Primary Tumor T Stage_T3a",
    "Neoplasm American Joint Committee on Cancer Clinical Primary Tumor T Stage_T3b",
    "Neoplasm American Joint Committee on Cancer Clinical Primary Tumor T Stage_T4",
    "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code_N1",
    "Positive Finding Lymph Node Hematoxylin and Eosin Staining Microscopy Count",
]


# ---------------------------------------------------------------------------
# 1. Stability selection (replaces PSO)
# ---------------------------------------------------------------------------
def stability_selection(
    X: pd.DataFrame,
    y: pd.Series,
    n_boot: int = 100,
    C: float = 0.1,
    threshold: float = 0.6,
    max_features: int = 15,
    random_state: int = RANDOM_STATE,
) -> Tuple[pd.Series, List[str]]:
    """Bootstrap-frequency L1 logistic selection.

    Runs L1-penalized logistic regression on `n_boot` stratified bootstrap
    resamples and records how often each feature gets a non-zero coefficient.
    Features selected in >= `threshold` of resamples are kept (capped at
    `max_features` by frequency).

    Returns
    -------
    (selection_frequency, selected_features)
    """
    rng = np.random.RandomState(random_state)
    freq = pd.Series(0.0, index=X.columns, dtype=float)

    for b in range(n_boot):
        idx = resample(np.arange(len(X)), stratify=y, random_state=rng.randint(1e9))
        Xb, yb = X.iloc[idx], y.iloc[idx]
        # guard: bootstrap sample must contain both classes
        if yb.nunique() < 2:
            continue
        # sklearn >=1.8: sparsity via l1_ratio=1 (liblinear), not penalty='l1'
        l1 = LogisticRegression(
            l1_ratio=1.0, solver="liblinear", C=C, max_iter=2000,
            class_weight="balanced", random_state=random_state,
        )
        l1.fit(Xb, yb)
        freq += (np.abs(l1.coef_[0]) > 1e-6).astype(float)

    freq /= max(n_boot, 1)
    ranked = freq.sort_values(ascending=False)
    selected = [c for c in ranked.index if ranked[c] >= threshold][:max_features]
    logger.info(
        "Stability selection: %d/%d features >= %.2f frequency (cap %d)",
        len(selected), len(X.columns), threshold, max_features,
    )
    return freq, selected


# ---------------------------------------------------------------------------
# 2. Calibrated Elastic Net (primary genomic model)
# ---------------------------------------------------------------------------
def build_elastic_net(random_state: int = RANDOM_STATE) -> CalibratedClassifierCV:
    """Calibrated Elastic Net logistic regression.

    Inner 5-fold CV picks C and l1_ratio by ROC-AUC on balanced weights;
    isotonic calibration on top restores interpretable probabilities so the
    0.5 threshold (rather than 0.0398) becomes usable.
    """
    base = LogisticRegressionCV(
        Cs=np.logspace(-3, 2, 15),
        l1_ratios=[0.1, 0.5, 0.9],
        solver="saga",
        use_legacy_attributes=False,
        cv=StratifiedKFold(5, shuffle=True, random_state=random_state),
        scoring="roc_auc",
        class_weight="balanced",
        max_iter=5000,
        random_state=random_state,
        n_jobs=-1,
    )
    return CalibratedClassifierCV(base, method="isotonic", cv=5)


# ---------------------------------------------------------------------------
# 3. Honest nested CV with selection inside each fold
# ---------------------------------------------------------------------------
@dataclass
class NestedCVResult:
    outer_aucs: List[float] = field(default_factory=list)
    outer_aps: List[float] = field(default_factory=list)
    fold_selected_features: List[List[str]] = field(default_factory=list)
    feature_stability: Optional[pd.Series] = None

    @property
    def mean_auc(self) -> float:
        return float(np.mean(self.outer_aucs))

    @property
    def std_auc(self) -> float:
        return float(np.std(self.outer_aucs))

    def summary(self) -> Dict[str, float]:
        return {
            "mean_outer_auc": self.mean_auc,
            "std_outer_auc": self.std_auc,
            "mean_outer_ap": float(np.mean(self.outer_aps)),
            "n_folds": len(self.outer_aucs),
        }


def nested_cv_with_selection(
    X: pd.DataFrame,
    y: pd.Series,
    candidate_features: Optional[Sequence[str]] = None,
    n_splits: int = 5,
    n_repeats: int = 1,
    n_boot: int = 50,
    stability_threshold: float = 0.6,
    max_features: int = 12,
    random_state: int = RANDOM_STATE,
) -> NestedCVResult:
    """Nested CV where stability selection happens inside every outer fold.

    This is the correction for the current pipeline's optimistic CV estimate:
    there, features were selected once on the full training data, so CV folds
    had already 'seen' the selection. Here each outer test fold is scored by a
    model whose features were chosen only from the outer-train partition.

    Parameters
    ----------
    X : preprocessed feature matrix (columns already filtered to candidates)
    candidate_features : restrict to this pool (e.g. MI top-200 ∩ external platform)
    """
    if candidate_features is not None:
        X = X[list(candidate_features)]

    cv = RepeatedStratifiedKFold(
        n_splits=n_splits, n_repeats=n_repeats, random_state=random_state
    )
    result = NestedCVResult()
    all_selected: List[str] = []

    for fold_i, (tr, te) in enumerate(cv.split(X, y)):
        X_tr, y_tr = X.iloc[tr], y.iloc[tr]
        X_te, y_te = X.iloc[te], y.iloc[te]

        _, feats = stability_selection(
            X_tr, y_tr, n_boot=n_boot, threshold=stability_threshold,
            max_features=max_features, random_state=random_state + fold_i,
        )
        if len(feats) == 0:  # degenerate fold: fall back to top-MI single gene
            feats = [X_tr.var().idxmax()]

        model = build_elastic_net(random_state + fold_i)
        model.fit(X_tr[feats], y_tr)
        prob = model.predict_proba(X_te[feats])[:, 1]

        result.outer_aucs.append(roc_auc_score(y_te, prob))
        result.outer_aps.append(average_precision_score(y_te, prob))
        result.fold_selected_features.append(feats)
        all_selected.extend(feats)

        logger.info(
            "Fold %2d: %2d features | AUC %.3f | AP %.3f",
            fold_i, len(feats), result.outer_aucs[-1], result.outer_aps[-1],
        )

    result.feature_stability = (
        pd.Series(all_selected).value_counts(normalize=True).sort_values(ascending=False)
    )
    return result


# ---------------------------------------------------------------------------
# 4. Clinical-only baseline (the reference to beat)
# ---------------------------------------------------------------------------
def clinical_baseline_columns(X: pd.DataFrame) -> List[str]:
    """Return clinical baseline columns actually present in X."""
    return [c for c in CLINICAL_BASELINE_FEATURES if c in X.columns]


def evaluate_clinical_baseline(
    X: pd.DataFrame,
    y: pd.Series,
    random_state: int = RANDOM_STATE,
) -> Dict[str, float]:
    """Repeated-CV AUC of the clinical-only model.

    A genomic model is only publishable if it beats this number (or matches it
    with added utility such as fewer inputs / better calibration).
    """
    cols = clinical_baseline_columns(X)
    if len(cols) < 3:
        logger.warning("Clinical baseline: only %d columns found", len(cols))
    model = build_elastic_net(random_state)
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=4, random_state=random_state)
    aucs, aps = [], []
    Xc = X[cols]
    for tr, te in cv.split(Xc, y):
        model.fit(Xc.iloc[tr], y.iloc[tr])
        prob = model.predict_proba(Xc.iloc[te])[:, 1]
        aucs.append(roc_auc_score(y.iloc[te], prob))
        aps.append(average_precision_score(y.iloc[te], prob))
    return {
        "clinical_auc_mean": float(np.mean(aucs)),
        "clinical_auc_std": float(np.std(aucs)),
        "clinical_ap_mean": float(np.mean(aps)),
        "n_features": len(cols),
    }


# ---------------------------------------------------------------------------
# 5. Cross-platform external validation helpers
# ---------------------------------------------------------------------------
def patient_zscore(df: pd.DataFrame, genes: Sequence[str]) -> pd.DataFrame:
    """Per-patient z-score across the given genes.

    Row-wise standardization removes per-array / per-library scale effects
    between RNA-Seq and microarray, leaving only relative expression --
    the part that can transfer between platforms.
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
    """Location-scale (ComBat-style) harmonization with frozen parameters.

    Estimates per-gene location/scale on `source` (training platform),
    standardizes `target` with its OWN batch statistics, then re-applies the
    source parameters -- the standard frozen ComBat transfer. No target
    labels are touched, so validation stays honest.
    """
    src = source[list(genes)].astype(float)
    tgt = target[list(genes)].astype(float)

    src_mu, src_sd = src.mean(axis=0), src.std(axis=0, ddof=1).replace(0, 1e-9)
    tgt_mu, tgt_sd = tgt.mean(axis=0), tgt.std(axis=0, ddof=1).replace(0, 1e-9)

    standardized = tgt.sub(tgt_mu, axis=1).div(tgt_sd, axis=1)
    return standardized.mul(src_sd, axis=1).add(src_mu, axis=1)


def common_gene_space(
    train_cols: Sequence[str],
    external_cols: Sequence[str],
) -> List[str]:
    """Intersection of gene columns across platforms (plain symbol match)."""
    external_set = set(external_cols)
    return [c for c in train_cols if c in external_set]


def recalibrate_threshold_for_prevalence(
    threshold: float, train_prevalence: float, external_prevalence: float
) -> float:
    """Shift a probability threshold for a cohort with different prevalence.

    Applies the logit-offset correction logit(p_ext) = logit(p_train) +
    logit(prev_ext) - logit(prev_train), valid under recalibrated sampling.
    """
    eps = 1e-6
    t = float(np.clip(threshold, eps, 1 - eps))
    logit_t = np.log(t / (1 - t))
    offset = (
        np.log(external_prevalence / (1 - external_prevalence))
        - np.log(train_prevalence / (1 - train_prevalence))
    )
    p = 1 / (1 + np.exp(-(logit_t + offset)))
    return float(np.clip(p, eps, 1 - eps))


# ---------------------------------------------------------------------------
# 6. End-to-end external validation on the common feature space
# ---------------------------------------------------------------------------
def cross_platform_validation(
    X_train: pd.DataFrame,
    y_train: pd.Series,
    X_ext: pd.DataFrame,
    y_ext: pd.Series,
    n_boot: int = 50,
    stability_threshold: float = 0.6,
    max_features: int = 12,
    normalization: str = "patient_zscore",  # or "frozen_combat"
    random_state: int = RANDOM_STATE,
) -> Dict[str, object]:
    """Train on the common gene space, validate on the external cohort.

    Fixes the two fatal flaws of notebooks 08/09:
      * selection now happens only among genes present on BOTH platforms;
      * row-wise / frozen-batch normalization removes platform scale before
        the model transfers.
    """
    genes = common_gene_space(X_train.columns, X_ext.columns)
    logger.info("Common gene space: %d genes", len(genes))
    if len(genes) < 5:
        raise ValueError(
            f"Only {len(genes)} common genes -- external validation not meaningful"
        )

    if normalization == "patient_zscore":
        Xtr = patient_zscore(X_train, genes)
        Xex = patient_zscore(X_ext, genes)
    elif normalization == "frozen_combat":
        Xex = frozen_combat(X_train, X_ext, genes)
        Xtr = X_train[genes].astype(float)
    else:  # raw columns
        Xtr, Xex = X_train[genes].astype(float), X_ext[genes].astype(float)

    # selection restricted to the transferable space, inside CV
    _, feats = stability_selection(
        Xtr, y_train, n_boot=n_boot, threshold=stability_threshold,
        max_features=max_features, random_state=random_state,
    )
    if not feats:
        feats = list(Xtr.std().sort_values(ascending=False).index[:max_features])

    model = build_elastic_net(random_state)
    model.fit(Xtr[feats], y_train)
    prob_ext = model.predict_proba(Xex[feats])[:, 1]

    train_prev = float(y_train.mean())
    ext_prev = float(y_ext.mean())
    thr05 = 0.5
    thr_adj = recalibrate_threshold_for_prevalence(thr05, train_prev, ext_prev)

    y_pred = (prob_ext >= thr_adj).astype(int)
    tp = int(((y_pred == 1) & (y_ext.values == 1)).sum())
    fn = int(((y_pred == 0) & (y_ext.values == 1)).sum())
    tn = int(((y_pred == 0) & (y_ext.values == 0)).sum())
    fp = int(((y_pred == 1) & (y_ext.values == 0)).sum())

    return {
        "n_common_genes": len(genes),
        "selected_features": feats,
        "external_auc": roc_auc_score(y_ext, prob_ext),
        "external_ap": average_precision_score(y_ext, prob_ext),
        "threshold_raw": thr05,
        "threshold_prevalence_adjusted": thr_adj,
        "sensitivity": tp / max(tp + fn, 1),
        "specificity": tn / max(tn + fp, 1),
        "tp": tp, "fn": fn, "tn": tn, "fp": fp,
        "normalization": normalization,
    }
