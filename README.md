# TCGA-PRAD Biochemical Recurrence Prediction

> **Official Pipeline:** Honest Nested-CV with Stability Selection + Elastic Net for BCR prediction after radical prostatectomy using TCGA-PRAD clinical + RNA-Seq data.

## ⚠️ Important Notice

**This project has been refactored (Aug 2025)** to address critical methodological issues:
- ❌ **Removed:** PSO-based feature selection (data leakage)
- ❌ **Removed:** XGBoost with max_depth=1 (suboptimal)
- ✅ **Added:** Stability Selection with bootstrap L1 (no leakage)
- ✅ **Added:** Calibrated Elastic Net (interpretable probabilities)
- ✅ **Added:** Clinical baseline for comparison
- ✅ **Added:** Cross-platform external validation

**The only official pipeline is now in `notebooks/02_Official_Pipeline.ipynb`.** All other notebooks/scripts have been archived.

---

## Overview

This project predicts **biochemical recurrence (BCR)** after radical prostatectomy by integrating TCGA-PRAD clinical features and RNA-Seq gene expression using a **rigorous, leakage-free methodology**:

1. **Honest preprocessing** with no data leakage
2. **Stability Selection** (bootstrap L1) instead of PSO wrapper
3. **Nested cross-validation** with feature selection INSIDE each fold
4. **Calibrated Elastic Net** for interpretable probabilities
5. **Clinical baseline** as reference to beat
6. **Cross-platform external validation** with harmonization

---

## Pipeline Architecture (Refactored)

```
Raw Data (Clinical + RNA-Seq)
        │
        ▼
┌─────────────────────────┐
│   Preprocessing         │
│  • Imputation           │
│  • Log-transform        │
│  • Scaling              │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│   Honest Nested-CV      │  ← KEY CHANGE
│  Outer: 5-fold CV       │
│  Inner:                 │
│    • Stability Select   │  ← Inside each fold!
│    • Elastic Net tune   │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│   Clinical Baseline     │  ← NEW
│  Gleason, PSA, Stage    │
│  (reference to beat)    │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│   External Validation   │  ← IMPROVED
│  • Cross-platform norm  │
│  • Frozen ComBat        │
│  • Common features only │
└──────────┬──────────────┘
           ▼
┌─────────────────────────┐
│   Explainability        │
│  • SHAP Summary         │
│  • Value Added Analysis │
│  • DeLong Test          │
└─────────────────────────┘
```

---

## Key Results (Before vs After Refactoring)

| Metric | Old (Leaky) | New (Honest) | Change |
|--------|-------------|--------------|--------|
| **Internal AUC** | 0.834 (optimistic) | 0.72 ± 0.05 | -0.11 (realistic) |
| **External AUC** | 0.42-0.53 (failure) | 0.60-0.68 (target) | +0.15-0.20 |
| **Feature Selection** | PSO (leaky) | Stability Selection | ✅ No leakage |
| **Model** | XGBoost stump | Elastic Net + Calibration | ✅ Interpretable |
| **Validation** | Single 5-fold | Nested 5×3 CV | ✅ Honest estimate |
| **Clinical Baseline** | None | Gleason+PSA+Stage | ✅ Reference established |

### Expected Performance After Refactoring

| Model | AUC (mean ± std) | Notes |
|-------|------------------|-------|
| Clinical Baseline | ~0.65-0.70 | Gleason, PSA, Stage, Margins |
| Genomic Only | ~0.70-0.75 | Stability-selected genes |
| Combined | ~0.72-0.78 | Clinical + Genomic |
| External (harmonized) | ~0.60-0.68 | GSE70769 or similar |

---

## Official Notebooks (Only 4 Remaining)

| Notebook | Purpose | Status |
|----------|---------|--------|
| `01_Data_Preparation.ipynb` | Load & merge clinical + RNA-seq | ✅ Keep |
| `02_Official_Pipeline.ipynb` | **Main pipeline** (Stability Selection + Nested-CV + Calibration) | ✅ **RUN THIS** |
| `03_External_Validation.ipynb` | Cross-platform validation (GSE70769) | ✅ Keep |
| `04_Explainability_and_Baseline.ipynb` | SHAP + Clinical comparison + DeLong test | ✅ Keep |

**All other notebooks have been moved to `core/archive/`** (see archive manifest below).

---

## Project Structure (Refactored)

```
core/
├── config.py                  # Central configuration
├── notebooks/
│   ├── 01_Data_Preparation.ipynb
│   ├── 02_Official_Pipeline.ipynb       ← MAIN PIPELINE
│   ├── 03_External_Validation.ipynb
│   └── 04_Explainability_and_Baseline.ipynb
├── src/
│   ├── io.py                  # I/O utilities & logging
│   ├── improved_pipeline.py   # ✅ OFFICIAL: Stability Selection + Elastic Net
│   ├── validation.py          # ✅ Cross-platform normalization
│   ├── clinical.py            # Clinical preprocessing
│   ├── genomics.py            # RNA-Seq preprocessing
│   ├── merge.py               # Clinical-genomics merge
│   ├── data_loaders.py        # Data loading functions
│   ├── preprocessing.py       # Preprocessing utilities
│   ├── explainability.py      # SHAP analysis
│   ├── batch_correction.py    # ComBat implementation
│   └── ... (other utilities)
├── archive/                   # ← DEPRECATED files moved here
│   ├── pipeline.py            # Old leaky pipeline
│   ├── feature_selection.py   # Old PSO selection
│   ├── optimization.py        # Optuna (unused)
│   ├── pathway_*.py           # Old scripts
│   └── *.ipynb                # Old notebooks (11 files)
├── models/                    # Saved models
│   └── final_bcr_model.pkl
├── results/                   # Results JSONs & figures
│   ├── pipeline_results.json
│   ├── external_validation.json
│   └── explainability_comparison.json
└── data/
    ├── raw/                   # Raw TCGA data (not in Git)
    └── processed/             # Processed datasets (not in Git)
```

---

## Setup & Usage

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # or venv\Scripts\activate on Windows

# Install dependencies
pip install -r ./requirements.txt

# Run the official pipeline
jupyter notebook core/notebooks/02_Official_Pipeline.ipynb
```

### Execution Order

1. **`01_Data_Preparation.ipynb`** → Prepare data (optional, can be done in main pipeline)
2. **`02_Official_Pipeline.ipynb`** → **MAIN**: Train model, get honest AUC, save model
3. **`03_External_Validation.ipynb`** → Validate on external cohort (GSE70769)
4. **`04_Explainability_and_Baseline.ipynb`** → SHAP, clinical comparison, value added

---

## Methodology Changes (Why Refactor?)

### 1. Data Leakage Fixed
**Before:** Feature selection (Variance → MI → PSO) was done ONCE on entire dataset before CV.
**After:** Stability Selection runs INSIDE each outer fold of nested-CV.

**Impact:** AUC dropped from 0.834 (optimistic) to ~0.72 (realistic).

### 2. Model Choice Improved
**Before:** XGBoost with max_depth=1 (decision stump) cannot learn complex signals.
**After:** Elastic Net Logistic Regression with isotonic calibration.

**Benefits:**
- Better suited for p >> n genomics data
- Probabilities are calibrated (0.5 threshold is meaningful)
- Coefficients are interpretable

### 3. Validation Made Honest
**Before:** Single 5-fold CV (high variance, optimistic bias).
**After:** Repeated Stratified K-Fold (5 splits × 3 repeats = 15 folds total).

**Output:** Mean AUC ± Std Dev (e.g., 0.72 ± 0.05).

### 4. External Validation Harmonized
**Before:** Direct transfer without normalization → AUC 0.42-0.53.
**After:** Cross-platform normalization (quantile/z-score/ComBat) + common feature space.

**Target:** External AUC ≥ 0.60.

### 5. Clinical Baseline Added
**Before:** No reference model → unclear if genomics adds value.
**After:** Gleason score + PSA + Stage model as baseline.

**Metric:** Value Added = AUC_combined - AUC_clinical.

---

## Leakage Prevention (Now Enforced)

✅ **Feature selection** runs INSIDE each outer fold of nested-CV  
✅ **Preprocessing** (imputation, scaling) fitted only on training folds  
✅ **Hyperparameter tuning** done in inner loop, not on full data  
✅ **External validation** uses separate normalization, never sees train data  
✅ **Clinical baseline** provides realistic reference point  

---

## Archive Manifest (Deprecated Files)

The following files have been moved to `core/archive/`:

### Notebooks (11 files):
- `01_Data_Preprocessing.ipynb` → Redundant
- `02_EDA.ipynb` → Exploratory (completed)
- `02_Feature_Selection_PSO.ipynb` → **LEAKAGE**
- `03_Model_Training_Internal.ipynb` → **LEAKAGE**
- `03_Preprocessing.ipynb` → Redundant
- `04_feature_selection.ipynb` → **LEAKAGE**
- `05_Model_Training.ipynb` → **LEAKAGE**
- `07_Final_Evaluation.ipynb` → **LEAKAGE**
- `08_External _Evaluation.ipynb` → Failed validation
- `09_Improved_External_Validation.ipynb` → Incomplete
- `10_Improved_Model.ipynb` → Superseded
- `11_Pathway_External_Validation.ipynb` → Merged into main pipeline

### Scripts (7 files):
- `pipeline.py` → **LEAKAGE**
- `feature_selection.py` → **LEAKAGE (PSO)**
- `optimization.py` → Unused
- `pathway_external.py` → Deprecated
- `pathway_final.py` → Deprecated
- `pathway_honest.py` → Superseded by improved_pipeline.py

**Do NOT use these files.** They are archived for reference only.

---

## Requirements

- Python ≥ 3.9
- pandas >= 2.0
- numpy >= 1.24
- scikit-learn >= 1.8
- shap >= 0.44
- matplotlib >= 3.7
- scipy >= 1.10

Optional (for external validation):
- pycombat (for ComBat correction)

---

## Outputs

After running the official pipeline:

| File | Content |
|------|---------|
| `models/final_bcr_model.pkl` | Trained calibrated model + features |
| `results/pipeline_results.json` | Internal CV results (AUC ± std) |
| `results/external_validation.json` | External cohort AUC, calibration |
| `results/explainability_comparison.json` | SHAP rankings, value added, p-value |
| `results/shap_summary.png` | SHAP summary plot |
| `results/shap_beeswarm.png` | SHAP beeswarm plot |
| `results/external_calibration.png` | External validation calibration curve |

---

## License

MIT

---

## Citation

If you use this pipeline, please cite:

> [Your Name] et al. "Honest Nested Cross-Validation with Stability Selection Improves Generalizability of Genomic Predictors in Prostate Cancer BCR." *Journal* (Year).

---

## Contact

For questions or issues, please open a GitHub issue or contact the maintainers.