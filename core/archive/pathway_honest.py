"""Honest pathway selection: rank by INTERNAL train CV only, evaluate external once."""
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

def pscores(df, sets, common):
    s = pd.DataFrame(index=df.index)
    for name, genes in sets.items():
        p = [g for g in genes if g in common]
        if len(p) >= 2:
            s[name] = df[p].mean(axis=1)
    return s

ps_t = pscores(X_train, ALL_PATHWAYS, wc)
ps_e = pscores(X_ext, ALL_PATHWAYS, wc)
ps_test = pscores(X_test, ALL_PATHWAYS, wc)

# ============================================================
# HONEST SELECTION: rank pathways by INTERNAL train CV-AUC only
# ============================================================
print('Internal (train 5-fold CV) AUC per pathway:')
internal_rank = {}
for col in ps_t.columns:
    lr1 = LogisticRegression(C=1, max_iter=2000, class_weight='balanced')
    cv_auc = cross_val_score(lr1, ps_t[[col]], y_train, cv=5, scoring='roc_auc').mean()
    internal_rank[col] = cv_auc
    print(f'  {col:22s} {cv_auc:.3f}')

ranked_internal = [k for k, v in sorted(internal_rank.items(), key=lambda x: -x[1])]
print(f'\nInternal ranking: {ranked_internal}')

# External evaluation ONCE per honest config
results = []
for k in [3, 5, 7, 10, 15]:
    cols = ranked_internal[:k]
    # tune C internally only
    best_C, best_cv = None, -1
    for C in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1]:
        lr = LogisticRegression(C=C, max_iter=2000, class_weight='balanced')
        cv = cross_val_score(lr, ps_t[cols], y_train, cv=5, scoring='roc_auc').mean()
        if cv > best_cv:
            best_cv, best_C = cv, C
    lr_final = LogisticRegression(C=best_C, max_iter=2000, class_weight='balanced')
    lr_final.fit(ps_t[cols], y_train)
    ae = roc_auc_score(y_ext, lr_final.predict_proba(ps_e[cols])[:, 1])
    at = roc_auc_score(y_test, lr_final.predict_proba(ps_test[cols])[:, 1])
    results.append((f'Honest Top{k} C={best_C}', ae, at, best_cv, k))
    print(f'Honest Top{k} (C={best_C}): internal CV={best_cv:.3f} | test={at:.3f} | external={ae:.3f}')

# ============================================================
# Bootstrap CI for the honest best
# ============================================================
best = max(results, key=lambda x: x[3])  # choose by INTERNAL CV, not external!
print(f'\nSelected by internal CV: {best[0]}')
k = best[4]
cols = ranked_internal[:k]
C_best = float(best[0].split('C=')[1])

lr_final = LogisticRegression(C=C_best, max_iter=2000, class_weight='balanced')
lr_final.fit(ps_t[cols], y_train)
prob_ext = lr_final.predict_proba(ps_e[cols])[:, 1]
auc_ext = roc_auc_score(y_ext, prob_ext)

rng = np.random.RandomState(42)
n = len(y_ext)
boot = []
for b in range(3000):
    idx = rng.choice(n, size=n, replace=True)
    yb = y_ext.values[idx]
    if len(np.unique(yb)) < 2:
        continue
    boot.append(roc_auc_score(yb, prob_ext[idx]))
ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])
p_value = 2 * min(np.mean(np.array(boot) <= 0.5), np.mean(np.array(boot) >= 0.5))

print(f'\n{"="*60}')
print(f'HONEST PATHWAY MODEL -> GSE70769')
print(f'  Pathways ({k}): {cols}')
print(f'  C={C_best}')
print(f'  External AUC = {auc_ext:.3f}')
print(f'  95% CI (3000 bootstraps) = [{ci_lo:.3f}, {ci_hi:.3f}]')
print(f'  CI above 0.5: {"YES" if ci_lo > 0.5 else "NO"}')
print(f'  p-value vs AUC=0.5: {p_value:.2e}')
print(f'  Internal test AUC = {best[2]:.3f}')
print(f'  Internal CV AUC   = {best[3]:.3f}')
print(f'{"="*60}')

# Coefficients
print('Coefficients:')
for col, coef in sorted(zip(cols, lr_final.coef_[0]), key=lambda x: -abs(x[1])):
    print(f'  {col:22s} {coef:+.4f}')

# DeLong-style sensitivity: also report AllPath honest
pd.DataFrame(results, columns=['config', 'external_auc', 'test_auc', 'cv_auc', 'k']).to_csv(
    'outputs/tables/pathway_honest_results.csv', index=False)
print('\nSaved -> outputs/tables/pathway_honest_results.csv')
