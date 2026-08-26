"""
Pathway gene sets — single source of truth for the Prostate BCR project.

Every pathway / gene-set used anywhere (notebooks, scripts, evaluation)
is defined here.  Import from this module only; do NOT re-define gene lists
in notebooks or scripts.

References:
  - Prolaris (Prolaris, Myriad Genetics): 46-gene cell-cycle progression score
  - Decipher (GenomeDx): 22-gene metastasis classifier
  - Androgen response: hallmark MSigDB HALLMARK_ANDROGEN_RESPONSE
  - EMT: hallmark MSigDB HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
  - Proliferation: MSigDB M5949 (proliferation markers)
"""

from __future__ import annotations

from typing import Dict, List

# ──────────────────────────────────────────────────────────────────────
# Individual pathway gene lists
# ──────────────────────────────────────────────────────────────────────

DECIPHER: List[str] = [
    "CEACAM1", "FLNA", "HES6", "KPNA2", "LCP1", "PLA2G7", "PTGER4",
    "RAB25", "SAA1", "SORD", "STOM", "TPX2", "TUBE1", "PDSS2",
    "SELENBP1", "SRD5A2", "TP53BP1",
]

PROLARIS: List[str] = [
    "BIRC5", "CDC20", "CDKN1A", "CENPF", "DUSP6", "EZH2", "FOXM1",
    "GTSE1", "KLK2", "KIF11", "KIF14", "KIF20A", "MCM2", "MCM5",
    "MCM7", "MKI67", "NDC80", "PCNA", "PLK1", "PTTG1", "RRM2",
    "SPP1", "TOP2A", "AURKA", "AURKB", "BUB1", "BUB1B", "CCNB1",
    "CCNB2", "CDCA3", "CDKN3", "CENPE", "CENPN", "DLGAP5", "EXO1",
    "GAS6", "HMMR", "KIF2C", "KIF4A", "MELK", "NCAPD2", "NUF2",
    "PBK", "RACGAP1", "RFC4", "TK1", "UBE2C", "ZWINT",
]

AR_SIGNALING: List[str] = [
    "AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "STEAP2", "ACPP", "CAMKK2",
]

EMT: List[str] = [
    "CDH1", "VIM", "CDH2", "SNAI1", "SNAI2", "ZEB1", "FN1", "CD44", "ITGA6",
]

PROLIFERATION: List[str] = [
    "MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "MCM7", "AURKA", "BIRC5", "CCNB1",
]

DNA_REPAIR: List[str] = [
    "BRCA1", "BRCA2", "RAD51", "ATM", "CHEK2", "XRCC2", "PARP1",
    "PALB2", "RAD54L", "GEN1",
]

PI3K_AKT: List[str] = [
    "PTEN", "PIK3CA", "AKT1", "MTOR", "RPS6KB1", "EIF4EBP1",
    "PDK1", "TSC1", "TSC2",
]

ANDROGEN_RESPONSE: List[str] = [
    "AR", "KLK3", "KLK2", "TMPRSS2", "SRD5A2", "HSD3B1", "CYP17A1",
]

CELL_CYCLE: List[str] = [
    "CCND1", "CCNE1", "CDK2", "CDK4", "CDK6", "RB1", "E2F1",
    "TP53", "CDKN1A", "CDKN2A", "CDKN1B",
]

STROMA: List[str] = [
    "ACTA2", "COL1A1", "COL3A1", "FAP", "PDGFRB", "POSTN",
    "THBS1", "TGFBI", "TAGLN", "VIM",
]

IMMUNE: List[str] = [
    "CD68", "CD8A", "CD4", "FOXP3", "PDCD1", "CTLA4",
    "LAG3", "CD274", "IFNG", "GZMB",
]

HYPOXIA: List[str] = [
    "HIF1A", "VEGFA", "CA9", "EGLN1", "EGLN3", "SLC2A1",
    "LOX", "P4HA1", "LDHA",
]

METABOLISM: List[str] = [
    "SLC2A1", "HK2", "PKM", "LDHA", "ACLY", "FASN",
    "SCD", "ACACA", "HMGCS2", "CPT1A",
]

WNT_BETA_CATENIN: List[str] = [
    "CTNNB1", "APC", "AXIN2", "LEF1", "TCF7", "MYC",
    "CCND1", "DVL2", "FRAT1", "WNT5A",
]

STRESS_RESPONSE: List[str] = [
    "HSP90AA1", "HSP90AB1", "HSPA5", "HSPA8", "HSPB1",
    "HSPD1", "HSPE1", "DNAJA1", "DNAJB1", "STIP1",
]

# PSA pathway (from features_config.py)
PSA_PATHWAY: List[str] = [
    "KLK3", "KLK2", "ACPP", "TMPRSS2", "AR", "NKX3-1", "STEAP2",
]


# ──────────────────────────────────────────────────────────────────────
# Master registry
# ──────────────────────────────────────────────────────────────────────

ALL_PATHWAYS: Dict[str, List[str]] = {
    "Decipher": DECIPHER,
    "Prolaris": PROLARIS,
    "AR_Signaling": AR_SIGNALING,
    "EMT": EMT,
    "Proliferation": PROLIFERATION,
    "DNA_Repair": DNA_REPAIR,
    "PI3K_AKT": PI3K_AKT,
    "Androgen_Response": ANDROGEN_RESPONSE,
    "Cell_Cycle": CELL_CYCLE,
    "Stroma": STROMA,
    "Immune": IMMUNE,
    "Hypoxia": HYPOXIA,
    "Metabolism": METABOLISM,
    "WNT_Beta_Catenin": WNT_BETA_CATENIN,
    "Stress_Response": STRESS_RESPONSE,
}

# Literature-motivated BCR-signature subset (pre-specified)
# These are published BCR-prognostic programs: Prolaris CCP, Decipher,
# cell cycle, DNA repair, proliferation.
BCR_LITERATURE_PATHWAYS: List[str] = [
    "Prolaris", "Decipher", "Cell_Cycle", "DNA_Repair", "Proliferation",
]

MIN_GENES_FOR_PATHWAY: int = 3


# ──────────────────────────────────────────────────────────────────────
# Pathway score computation
# ──────────────────────────────────────────────────────────────────────

def compute_pathway_scores(
    df,
    pathway_dict: Dict[str, List[str]] | None = None,
    common_genes: set | None = None,
    min_genes: int = MIN_GENES_FOR_PATHWAY,
):
    """Compute mean-expression pathway scores for a gene expression DataFrame.

    Parameters
    ----------
    df : DataFrame with genes as columns.
    pathway_dict : Mapping of pathway name → gene list. Defaults to ALL_PATHWAYS.
    common_genes : If provided, only use genes from this intersection.
    min_genes : Minimum genes required to compute a pathway score.

    Returns
    -------
    DataFrame with one column per pathway (mean of available genes).
    """
    import pandas as pd

    if pathway_dict is None:
        pathway_dict = ALL_PATHWAYS

    df_cols = set(df.columns)
    scores = pd.DataFrame(index=df.index)

    for name, genes in pathway_dict.items():
        available = [g for g in genes if g in df_cols]
        if common_genes is not None:
            available = [g for g in available if g in common_genes]
        if len(available) >= min_genes:
            scores[name] = df[available].mean(axis=1)

    return scores
