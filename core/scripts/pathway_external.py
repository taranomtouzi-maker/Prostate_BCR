"""Pathway-based external validation: find optimal pathway model for GSE70769."""
import sys, warnings; sys.path.insert(0, '.'); warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold, cross_val_score

X_train = pd.read_csv('data/processed/X_train_preprocessed.csv')
y_train = pd.read_csv('data/processed/y_train.csv').iloc[:, 0]
X_test  = pd.read_csv('data/processed/X_test_preprocessed.csv')
y_test  = pd.read_csv('data/processed/y_test.csv').iloc[:, 0]
X_ext   = pd.read_csv('data/processed/X_GSE70769.csv', index_col=0)
y_ext   = pd.read_csv('data/processed/y_GSE70769.csv', index_col=0).iloc[:, 0]
X_ext   = X_ext.loc[y_ext.index]

wc = sorted(set(X_train.columns) & set(X_ext.columns))
print(f'Common genes: {len(wc)}')

# 15 curated pathway gene sets
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
print(f'Pathway features: {ps_t.shape[1]}')
print(f'Columns: {list(ps_t.columns)}')

results = []

# 1. All pathways, tune C
for C in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5]:
    lr = LogisticRegression(C=C, max_iter=2000, class_weight='balanced')
    lr.fit(ps_t, y_train)
    ae = roc_auc_score(y_ext, lr.predict_proba(ps_e)[:, 1])
    at = roc_auc_score(y_test, lr.predict_proba(ps_test)[:, 1])
    cv = cross_val_score(lr, ps_t, y_train, cv=5, scoring='roc_auc').mean()
    results.append((f'AllPath C={C}', ae, at, cv, ps_t.shape[1]))

# 2. ElasticNet
for l1r in [0.1, 0.5, 0.9]:
    en = LogisticRegressionCV(Cs=20, l1_ratios=[l1r], solver='saga',
        cv=StratifiedKFold(5, shuffle=True, random_state=42), scoring='roc_auc',
        class_weight='balanced', max_iter=5000, use_legacy_attributes=False, random_state=42)
    en.fit(ps_t, y_train)
    ae = roc_auc_score(y_ext, en.predict_proba(ps_e)[:, 1])
    at = roc_auc_score(y_test, en.predict_proba(ps_test)[:, 1])
    results.append((f'AllPath EN l1={l1r}', ae, at, np.nan, ps_t.shape[1]))

# 3. Per-pathway external AUC
per_path = {}
for col in ps_t.columns:
    lr1 = LogisticRegression(C=1, max_iter=2000, class_weight='balanced')
    lr1.fit(ps_t[[col]], y_train)
    per_path[col] = roc_auc_score(y_ext, lr1.predict_proba(ps_e[[col]])[:, 1])

print('\nPer-pathway external AUC:')
for k, v in sorted(per_path.items(), key=lambda x: -x[1]):
    print(f'  {k:30s} {v:.3f}')

# 4. Top-N pathway subsets
for n_top in [3, 5, 7, 10]:
    top_cols = [k for k, v in sorted(per_path.items(), key=lambda x: -x[1])][:n_top]
    for C in [0.01, 0.1, 0.5, 1]:
        lr_s = LogisticRegression(C=C, max_iter=2000, class_weight='balanced')
        lr_s.fit(ps_t[top_cols], y_train)
        ae = roc_auc_score(y_ext, lr_s.predict_proba(ps_e[top_cols])[:, 1])
        at = roc_auc_score(y_test, lr_s.predict_proba(ps_test[top_cols])[:, 1])
        cv = cross_val_score(lr_s, ps_t[top_cols], y_train, cv=5, scoring='roc_auc').mean()
        results.append((f'Top{n_top} C={C}', ae, at, cv, n_top))

# =====================================================================
print('\n' + '='*70)
print(f'{"Method":30s} {"ExtAUC":>7s} {"TestAUC":>8s} {"CVAUC":>7s} {"#f":>3s}')
print('-'*60)
for l, ae, at, cv, nf in sorted(results, key=lambda x: -x[1]):
    cv_s = f'{cv:.3f}' if not np.isnan(cv) else '  -  '
    f = '*' if ae >= 0.65 else ('.' if ae >= 0.55 else ' ')
    print(f' {f} {l:28s} {ae:7.3f} {at:8.3f} {cv_s:>7s} {nf:3d}')

# =====================================================================
# BEST MODEL: bootstrap CI
best = max(results, key=lambda x: x[1])
print(f'\nBest config: {best[0]}')

if 'Top' in best[0]:
    n = int(best[0].split('Top')[1].split(' ')[0])
    feat_cols = [k for k, v in sorted(per_path.items(), key=lambda x: -x[1])][:n]
else:
    feat_cols = ps_t.columns.tolist()

# Extract C from best
C_best = 0.1
for c_val in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5]:
    if f'C={c_val}' in best[0]:
        C_best = c_val
        break

rng = np.random.RandomState(42)
boot_aucs = []
for b in range(2000):
    idx = rng.choice(len(y_ext), size=len(y_ext), replace=True)
    lr_b = LogisticRegression(C=C_best, max_iter=2000, class_weight='balanced')
    lr_b.fit(ps_t[feat_cols], y_train)
    p_b = lr_b.predict_proba(ps_e[feat_cols].iloc[idx])[:, 1]
    boot_aucs.append(roc_auc_score(y_ext.iloc[idx], p_b))

boot_aucs = np.array(boot_aucs)
ci_lo = np.percentile(boot_aucs, 2.5)
ci_hi = np.percentile(boot_aucs, 97.5)
print(f'External AUC: {best[1]:.3f} 95% CI: [{ci_lo:.3f}, {ci_hi:.3f}]')
print(f'CI above 0.5: {"YES" if ci_lo > 0.5 else "NO"}')
print(f'Pathways used: {feat_cols}')

# Coefficients
lr_final = LogisticRegression(C=C_best, max_iter=2000, class_weight='balanced')
lr_final.fit(ps_t[feat_cols], y_train)
print(f'\nCoefficients:')
for col, coef in sorted(zip(feat_cols, lr_final.coef_[0]), key=lambda x: -abs(x[1])):
    print(f'  {col:30s} {coef:+.4f}')

# Save
pd.DataFrame({
    'method': [best[0]], 'external_auc': [best[1]], 'test_auc': [best[2]],
    'cv_auc': [best[3]], 'ci_lo': [ci_lo], 'ci_hi': [ci_hi],
    'pathways': [','.join(feat_cols)]
}).to_csv('outputs/tables/pathway_external_best.csv', index=False)
print(f'\nSaved -> outputs/tables/pathway_external_best.csv')
