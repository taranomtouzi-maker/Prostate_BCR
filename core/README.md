# TCGA-PRAD Biochemical Recurrence Prediction

> Multi-omics machine learning pipeline for predicting biochemical recurrence (BCR)
> after radical prostatectomy using TCGA-PRAD clinical + RNA-Seq data, with PSO-based
> feature selection and pathway-based external validation.

## Overview

This project predicts **biochemical recurrence (BCR)** after radical prostatectomy
by integrating TCGA-PRAD clinical features and RNA-Seq gene expression (~18,900 genes).
The pipeline implements:

1. **Rigorous preprocessing** with leakage auditing
2. **Two-stage feature selection**: Mutual Information → Binary PSO
3. **Nested cross-validation** model comparison across 6 classifiers
4. **Dual external validation**: gene-based + pathway-based on GSE70769

## Pipeline Architecture

```
Raw Data (Clinical + RNA-Seq)
        │
        ▼
┌─────────────────────────┐
│  NB01: Data Preprocessing│
│  • Clinical cleaning     │
│  • RNA-Seq preprocessing │
│  • Merge & leakage audit │
│  • Stratified train/test │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│  NB02: Feature Selection │
│  • Variance Threshold    │
│  • Mutual Information    │
│  • Feature Engineering   │
│  • Binary PSO (30 feat)  │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│  NB03: Model Training    │
│  • Nested 5-Fold CV      │
│  • 6-classifier compare  │
│  • XGBoost tuning        │
│  • Held-out test eval    │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│  NB04: Evaluation        │
│  • Internal recap        │
│  • Gene-ext validation   │
│  • Pathway-ext validation│
│  • GSE54460 triangulate  │
│  • Unified results table │
└─────────────────────────┘
```

## Key Results

### Internal Validation (TCGA-PRAD, held-out test set, n=86)

| Metric | Value | 95% CI |
|--------|-------|--------|
| **ROC-AUC** | **~0.82** | [0.70, 0.92] |
| PR-AUC | — | — |
| Sensitivity | 83.3% | — |
| Specificity | — | — |
| Balanced Accuracy | ~76% | — |
| Selected Features | 37 | — |

### External Validation (GSE70769, n=93, microarray)

| Approach | AUC | 95% CI | Notes |
|----------|-----|--------|-------|
| Gene-based (XGBoost) | ~0.55–0.65 | — | Honest baseline; platform noise degrades performance |
| Pathway-based (all 15) | ~0.65–0.71 | — | Robust to platform differences |
| Pathway-based (BCR-lit 5) | ~0.60–0.68 | — | Literature-curated subset |

> **Why pathway scores outperform gene-level features externally:**
> Pathway scores aggregate signal across multiple genes within biologically
> meaningful programs (cell cycle, AR signaling, proliferation). This
> aggregation averages out platform-specific noise while preserving
> conserved biological signal that transfers across RNA-Seq ↔ microarray.
> Individual gene expression levels are highly platform-dependent.

### Top Biomarkers Identified

| Rank | Gene | Importance |
|------|------|------------|
| 1 | DYNLT1 | 0.074 |
| 2 | POU2AF1 | 0.066 |
| 3 | SOCS2 | 0.058 |
| 4 | CNTRL | 0.040 |
| 5 | HSD11B1L | 0.035 |

## Project Structure

```
core/
├── config.py                  # Central configuration (seeds, paths, hyperparams)
├── requirements.txt           # Pinned dependencies
├── notebooks/
│   ├── 01_Data_Preprocessing.ipynb     # Clinical + RNA-Seq loading & merging
│   ├── 02_Feature_Selection_PSO.ipynb  # MI → PSO feature selection
│   ├── 03_Model_Training_Internal.ipynb# Nested CV, tuning, test eval
│   └── 04_Final_Evaluation_External.ipynb # Unified internal + external
├── src/
│   ├── io.py                  # I/O utilities & logging
│   ├── clinical.py            # Clinical preprocessing
│   ├── genomics.py            # RNA-Seq preprocessing
│   ├── merge.py               # Clinical-genomics merge
│   ├── leakage.py             # Leakage audit
│   ├── pathways.py            # ★ Pathway gene sets (single source of truth)
│   ├── data_loaders.py        # ★ External cohort loaders + normalization
│   ├── feature_engineering.py # ★ Engineered features (clinical + pathway)
│   ├── feature_selection.py   # MI + Binary PSO
│   ├── models.py              # Model factories & registry
│   ├── pipeline.py            # Nested CV & evaluation
│   ├── evaluation.py          # Metrics, bootstrap CIs, DCA
│   ├── visualization.py       # Plotting utilities
│   ├── explainability.py      # SHAP & biomarker ranking
│   ├── preprocessing.py       # Sklearn Pipeline transformers
│   ├── validation.py          # Feature validation helpers
│   ├── improved_pipeline.py   # Stability selection + ElasticNet (reference)
│   └── clinical_utility.py    # DCA, clinical impact curves
├── scripts/
│   └── triangulate_54460.py   # GSE54460 triangulation (optional)
├── data/
│   ├── raw/                   # Raw TCGA data (not in Git)
│   ├── processed/             # Processed datasets (not in Git)
│   └── external/              # External validation data (not in Git)
└── outputs/                   # Models, figures, tables (not in Git)
```

## Setup

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # or venv\Scripts\activate on Windows

# Install dependencies
pip install -r requirements.txt

# Place raw data files
# TCGA: data/raw/data_clinical_patient.tsv, data/raw/data_mrna_seq_v2_rsem.txt
# GSE70769: data/external/X_GSE70769.csv, data/external/y_GSE70769.csv
```

## Execution Order

Run notebooks **in order** — each produces artifacts consumed by the next:

```bash
# From the core/ directory:
jupyter notebook notebooks/01_Data_Preprocessing.ipynb
jupyter notebook notebooks/02_Feature_Selection_PSO.ipynb
jupyter notebook notebooks/03_Model_Training_Internal.ipynb
jupyter notebook notebooks/04_Final_Evaluation_External.ipynb
```

Or non-interactively:
```bash
jupyter nbconvert --execute notebooks/01_Data_Preprocessing.ipynb
jupyter nbconvert --execute notebooks/02_Feature_Selection_PSO.ipynb
jupyter nbconvert --execute notebooks/03_Model_Training_Internal.ipynb
jupyter nbconvert --execute notebooks/04_Final_Evaluation_External.ipynb
```

## Reproducibility

- **Random seed:** All seeds centralized in `config.RANDOM_STATE = 42`
- **No leakage:** Feature selection fitted only on training data; test data never accessed during selection or tuning
- **Logging:** All modules use `src.io.logger` (structured, timestamped)
- **Deterministic:** XGBoost `tree_method="hist"` with fixed seed

## External Validation Design

### Why Two Strategies?

1. **Gene-based:** Uses the exact same model features on external data (after patient-zscore normalization). Reports the honest transferability of the gene signature.

2. **Pathway-based:** Computes 15 curated pathway scores (mean expression of biologically coherent gene sets) and trains a logistic regression. More robust because:
   - Pathway scores aggregate multi-gene signal, averaging out platform noise
   - Biological pathways are evolutionarily conserved between RNA-Seq and microarray
   - Mean-expression scores are less sensitive to normalization differences
   - Literature-validated pathways (Prolaris, Decipher, cell cycle) are known BCR prognostic programs

### Cross-Platform Normalization

| Method | When Used | Description |
|--------|-----------|-------------|
| Quantile (default) | Gene-level external | Match expression distributions across platforms |
| Patient z-score | Gene-level transfer | Row-wise standardization removes per-array scale |
| Frozen ComBat | Optional | Location-scale harmonization with source parameters |

## Leaked Results Prevention

| Check | Status |
|-------|--------|
| Feature selection on full train | ⚠️ Acceptable for final model (test never seen) |
| Nested CV for honest estimate | ✅ NB03 uses outer CV |
| External: no target leakage | ✅ Only common genes + pre-specified pathways |
| External: C tuned on internal CV only | ✅ NB04 tunes C on TCGA, evaluates on GSE70769 |

## License

MIT
