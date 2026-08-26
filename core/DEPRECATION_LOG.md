# Deprecation Log — Consolidation Sprint

> Date: 2026-08-26
> Reason: Consolidation from 11 notebooks + 4 scripts → 4 core notebooks
> All logic from deprecated files has been absorbed into `src/` modules or the 4 core notebooks.

---

## Notebooks to Archive/Delete

| File | Status | Justification |
|------|--------|---------------|
| `02_EDA.ipynb` | **DELETE** | Exploratory analysis; not part of reproducible pipeline. No artifacts generated. |
| `03_Preprocessing.ipynb` | **MERGED → 01** | Preprocessing logic consolidated into `01_Data_Preprocessing.ipynb`. |
| `04_feature_selection.ipynb` | **MERGED → 02** | Feature selection logic consolidated into `02_Feature_Selection_PSO.ipynb`. |
| `05_Model_Training.ipynb` | **MERGED → 03** | Model training logic consolidated into `03_Model_Training_Internal.ipynb`. |
| `06_Explainability.ipynb` | **DELETE** | Explainability is post-hoc; SHAP analysis can be added as optional section in `03` or `04`. Core pipeline doesn't depend on it. |
| `07_Final_Evaluation.ipynb` | **MERGED → 04** | Internal evaluation logic merged into `04_Final_Evaluation_External.ipynb`. |
| `08_External _Evaluation.ipynb` | **DELETE** | Had parsing errors, missing clinical columns, no cross-platform normalization. Replaced by `04_Final_Evaluation_External.ipynb`. |
| `09_Improved_External_Validation.ipynb` | **DELETE** | Incremental fix of notebook 08; still incomplete. Superseded by `04_Final_Evaluation_External.ipynb`. |
| `10_Improved_Model.ipynb` | **DELETE** | Experimental notebook; stability selection + ElasticNet logic moved to `src/improved_pipeline.py` (kept as utility). |
| `11_Pathway_External_Validation.ipynb` | **DELETE** | Pathway logic fully consolidated into `src/pathways.py` and `04_Final_Evaluation_External.ipynb`. |

**After cleanup: 4 notebooks remain:**
1. `01_Data_Preprocessing.ipynb`
2. `02_Feature_Selection_PSO.ipynb`
3. `03_Model_Training_Internal.ipynb`
4. `04_Final_Evaluation_External.ipynb`

---

## Scripts to Archive/Delete

| File | Status | Justification |
|------|--------|---------------|
| `scripts/pathway_external.py` | **DELETE** | Pathway gene sets + evaluation logic consolidated into `src/pathways.py`. |
| `scripts/pathway_final.py` | **DELETE** | Duplicate of `pathway_honest.py` with minor variations. Logic in `src/pathways.py`. |
| `scripts/pathway_honest.py` | **DELETE** | "Honest" pathway selection logic moved to `src/pathways.py::honest_pathway_evaluation()`. |
| `scripts/triangulate_54460.py` | **ARCHIVE** | Useful for GSE54460 triangulation; kept as `scripts/triangulate_54460.py` but not part of core pipeline. |

---

## Source Modules — Merge/Delete

| File | Status | Justification |
|------|--------|---------------|
| `src/batch_correction.py` | **DELETE** | Overlaps with `src/validation.py` and `src/improved_pipeline.py`. ComBat logic moved to `src/data_loaders.py`. |
| `src/survival_analysis.py` | **ARCHIVE** | Useful for KM curves / C-index but not needed for BCR binary classification pipeline. Can be re-added for time-to-event analysis. |
| `src/optimization.py` | **ARCHIVE** | Optuna optimization not used in final pipeline. Kept for future reference. |
| `src/features_config.py` | **MERGED → pathways.py** | Gene set definitions merged into `src/pathways.py`. |
| `src/preprocessing.py` | **KEEP** | Used by notebooks; provides sklearn-compatible pipelines. |
| `src/validation.py` | **TRIM** | Cross-platform normalization moved to `src/data_loaders.py`. Keep feature validation functions only. |

**After cleanup: 14 src modules remain:**
- `__init__.py`, `batch_correction.py` (removed), `clinical.py`, `clinical_utility.py`, `data_loaders.py` (new), `evaluation.py`, `explainability.py`, `feature_engineering.py` (new, extracted from feature_selection.py), `feature_selection.py`, `features_config.py` (deprecated), `genomics.py`, `improved_pipeline.py`, `io.py`, `leakage.py`, `merge.py`, `models.py`, `optimization.py` (archived), `pathways.py` (new), `pipeline.py`, `preprocessing.py`, `survival_analysis.py` (archived), `validation.py` (trimmed), `visualization.py`

---

## Documentation to Update/Delete

| File | Status | Justification |
|------|--------|---------------|
| `PIPELINE_CONSISTENCY_FIX.md` | **DELETE** | Troubleshooting doc for old pipeline; no longer relevant. |
| `IMPROVEMENT_PLAN.md` | **KEEP** | Historical reference for why changes were made. |
| `doc.md` | **DELETE** | Old documentation; replaced by updated README.md. |
| `CLEANUP_REPORT.md` | **REPLACE** | This file replaces old cleanup report. |
| `CODE_AUDIT_REPORT.md` | **REPLACE** | This file supersedes old audit. |

---

## Migration Summary

**What was preserved:**
- `config.py` — Enhanced with new constants, backward-compatible
- `src/io.py` — Unchanged
- `src/clinical.py` — Unchanged
- `src/genomics.py` — Unchanged
- `src/merge.py` — Unchanged
- `src/leakage.py` — Unchanged
- `src/models.py` — Unchanged
- `src/pipeline.py` — Unchanged
- `src/feature_selection.py` — Kept; feature engineering extracted
- `src/visualization.py` — Kept; may add DCA plots
- `src/explainability.py` — Kept
- `src/preprocessing.py` — Kept
- `src/evaluation.py` — Enhanced with DCA, bootstrap CI for all metrics
- `src/improved_pipeline.py` — Kept as reference; stability selection used in 04
- `src/clinical_utility.py` — Kept; DCA used in 04

**What was created:**
- `src/pathways.py` — Single source of truth for all pathway gene sets
- `src/data_loaders.py` — `load_gse70769()` + cross-platform normalization
- `src/feature_engineering.py` — Extracted from `feature_selection.py`
- 4 core notebooks (see above)
- `requirements.txt` — Pinned dependencies
- `DEPRECATION_LOG.md` — This file
