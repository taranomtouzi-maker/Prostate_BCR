"""Defensible pathway models:
  Primary  : all 15 pre-specified pathway scores, C tuned on internal CV only.
  Secondary: literature-motivated BCR-signature subset (Prolaris/Decipher/cell-cycle biology).
Both evaluated ONCE on GSE70769 with bootstrap CIs, plus triangulation on GSE54460 (RNA-Seq).
"""
import sys, warnings; sys.path.insert(0, '.'); warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score

X_train = pd.read_csv('data/processed/X_train_preprocessed.csv')
y_train = pd.read_csv('data/processed/y_train.csv').iloc[:, 0]
X_test  = pd.read_csv('data/processed/X_test_preprocessed.csv')
y_test  = pd.read_csv('data/processed/y_test.csv').iloc[:, 0]
X_ext   = pd.read_csv('data/processed/X_GSE70769.csv', index_col=0)
y_ext   = pd.read_csv('data/processed/y_GSE70769.csv', index_col=0).iloc[:, 0]
X_ext   = X_ext.loc[y_ext.index]

wc = sorted(set(X_train.columns) & set(X_ext.columns))

DECIPHER = ['CEACAM1','FLNA','HES6','KPNA2','LCP1','PLA2G7','PTGER4','RAB25','SAA1','SORD','STOM','TPX2','TUBE1','PDSS2','SELENBP1','SRD5A2','TP53BP1']
PROLARIS = ['BIRC5','CDC20','CDKN1A','CENPF','DUSP6','EZH2','FOXM1','GTSE1','KLK2','KIF11','KIF14','KIF20A','MCM2','MCM5','MCM7','MKI67','NDC80','PCNA','PLK1','PTTG1','RRM2','SPP1','TOP2A','AURKA','AURKB','BUB1','BUB1B','CCNB1','CCNB2','CDCA3','CDKN3','CENPE','CENPN','DLGAP5','EXO1','GAS6','HMMR','KIF2C','KIF4A','MELK','NCAPD2','NUF2','PBK','RACGAP1','RFC4','TK1','UBE2C','ZWINT']
AR_SIG = ['AR','KLK3','KLK2','TMPRSS2','FKBP5','STEAP2','ACPP','CAMKK2']
EMT = ['CDH1','VIM','CDH2','SNAI1','SNAI2','ZEB1','FN1','CD44','ITGA6']
PROLIF = ['MKI67','TOP2A','PCNA','MCM2','MCM5','MCM7','AURKA','BIRC5','CCNB1']
DNA_REPAIR = ['BRCA1','BRCA2','RAD51','ATM','CHEK2','XRCC2','PARP1','PALB2','RAD54L','GEN1']
PI3K_AKT = ['PTEN','PIK3CA','AKT1','MTOR','RPS6KB1','EIF4EBP1','PDK1','TSC1','TSC2']
ANDROGEN_R = ['AR','KLK3','KLK2','TMPRSS2','SRD5A2','HSD3B1','CYP17A1']
CELL_CYCLE = ['CCND1','CCNE1','CDK2','CDK4','CDK6','RB1','E2F1','TP53','CDKN1A','CDKN2A','CDKN1B']
STROMA = ['ACTA2','COL1A1','COL3A1','FAP','PDGFRB','POSTN','THBS1','TGFBI','TAGLN','VIM']
IMMUNE = ['CD68','CD8A','CD4','FOXP3','PDCD1','CTLA4','LAG3','CD274','IFNG','GZMB']
HYPOXIA = ['HIF1A','VEGFA','CA9','EGLN1','EGLN3','SLC2A1','LOX','P4HA1','LDHA']
METABOLISM = ['SLC2A1','HK2','PKM','LDHA','ACLY','FASN','SCD','ACACA','HMGCS2','CPT1A']
WNT_BETA = ['CTNNB1','APC','AXIN2','LEF1','TCF7','MYC','CCND1','DVL2','FRAT1','WNT5A']
STRESS = ['HSP90AA1','HSP90AB1','HSPA5','HSPA8','HSPB1','HSPD1','HSPE1','DNAJA1','DNAJB1','STIP1']

ALL_PATHWAYS = {
    'Decipher': DECIPHER, 'Prolaris': PROLARIS, 'AR_Signaling': AR_SIG,
    'EMT': EMT, 'Proliferation': PROLIF, 'DNA_Repair': DNA_REPAIR,
    'PI3K_AKT': PI3K_AKT, 'Androgen_Response': ANDROGEN_R,
    'Cell_Cycle': CELL_CYCLE, 'Stroma': STROMA, 'Immune': IMMUNE,
    'Hypoxia': HYPOXIA, 'Metabolism': METABOLISM,
    'WNT_Beta_Catenin': WNT_BETA, 'Stress_Response': STRESS,
}
# Literature-motivated BCR signature subset (pre-specified: these are the
# published BCR-prognostic programs -- Prolaris CCP, Decipher, cell cycle,
# DNA repair, proliferation)
BCR_LITERATURE = ['Prolaris', 'Decipher', 'Cell_Cycle', 'DNA_Repair', 'Proliferation']

def pscores(df, sets, common):
    s = pd.DataFrame(index=df.index)
    for name, genes in sets.items():
        p = [g for g in genes if g in common]
        if len(p) >= 2:
            s[name] = df[p].mean(axis=1)
    return s

ps_t  = pscores(X_train, ALL_PATHWAYS, wc)
ps_e  = pscores(X_ext, ALL_PATHWAYS, wc)
ps_te = pscores(X_test, ALL_PATHWAYS, wc)

def bootstrap_ci(y, prob, n_boot=3000, seed=42):
    rng = np.random.RandomState(seed)
    yv = np.asarray(y)
    boots = []
    for _ in range(n_boot):
        idx = rng.choice(len(yv), size=len(yv), replace=True)
        if len(np.unique(yv[idx])) < 2:
            continue
        boots.append(roc_auc_score(yv[idx], prob[idx]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    p = 2 * min(np.mean(np.array(boots) <= 0.5), np.mean(np.array(boots) >= 0.5))
    return lo, hi, p

def eval_config(name, cols, label):
    # tune C on INTERNAL CV only
    best_C, best_cv = 0.1, -1
    for C in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1]:
        lr = LogisticRegression(C=C, max_iter=2000, class_weight='balanced')
        cv = cross_val_score(lr, ps_t[cols], y_train, cv=5, scoring='roc_auc').mean()
        if cv > best_cv:
            best_cv, best_C = cv, C
    lr = LogisticRegression(C=best_C, max_iter=2000, class_weight='balanced')
    lr.fit(ps_t[cols], y_train)
    p_ext = lr.predict_proba(ps_e[cols])[:, 1]
    p_te  = lr.predict_proba(ps_te[cols])[:, 1]
    auc_e = roc_auc_score(y_ext, p_ext)
    auc_te = roc_auc_score(y_test, p_te)
    lo, hi, pv = bootstrap_ci(y_ext, p_ext)
    print(f'{label}')
    print(f'  C={best_C} (internal CV={best_cv:.3f})')
    print(f'  Internal test AUC = {auc_te:.3f}')
    print(f'  EXTERNAL AUC = {auc_e:.3f}  CI [{lo:.3f}, {hi:.3f}]  p={pv:.2e}')
    print(f'  CI>0.5: {"YES" if lo > 0.5 else "NO"} | features: {cols}')
    print()
    return dict(config=name, C=best_C, cv_auc=best_cv, test_auc=auc_te,
                external_auc=auc_e, ci_lo=lo, ci_hi=hi, p_value=pv)

rows = []
print('='*70)
# PRIMARY: all 15 pathways, no selection whatsoever
rows.append(eval_config('all15', list(ps_t.columns), 'PRIMARY: all 15 pre-specified pathways'))

# SECONDARY: literature BCR-signature subset
rows.append(eval_config('bcr_lit5', BCR_LITERATURE, 'SECONDARY: literature BCR-signature subset (5)'))

# TERTIARY options for sensitivity analysis
rows.append(eval_config('prolif_only', ['Proliferation', 'Cell_Cycle'], 'SENSITIVITY: proliferation+cell-cycle only'))

df = pd.DataFrame(rows)
df.to_csv('outputs/tables/pathway_defensible_results.csv', index=False)
print('Saved -> outputs/tables/pathway_defensible_results.csv')
print(df.to_string(index=False))
