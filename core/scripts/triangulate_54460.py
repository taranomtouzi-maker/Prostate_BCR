"""Triangulate the BCR-literature pathway model on GSE54460 (2nd external cohort, RNA-Seq FFPE)."""
import sys, warnings; sys.path.insert(0, '.'); warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score

X_train = pd.read_csv('data/processed/X_train_preprocessed.csv')
y_train = pd.read_csv('data/processed/y_train.csv').iloc[:, 0]
X54460 = pd.read_csv('data/interim/X_GSE54460.csv', index_col=0)  # log2(FPKM+1)
y54460 = pd.read_csv('data/interim/y_GSE54460.csv', index_col=0).iloc[:, 0]

DECIPHER = ['CEACAM1','FLNA','HES6','KPNA2','LCP1','PLA2G7','PTGER4','RAB25','SAA1','SORD','STOM','TPX2','TUBE1','PDSS2','SELENBP1','SRD5A2','TP53BP1']
PROLARIS = ['BIRC5','CDC20','CDKN1A','CENPF','DUSP6','EZH2','FOXM1','GTSE1','KLK2','KIF11','KIF14','KIF20A','MCM2','MCM5','MCM7','MKI67','NDC80','PCNA','PLK1','PTTG1','RRM2','SPP1','TOP2A','AURKA','AURKB','BUB1','BUB1B','CCNB1','CCNB2','CDCA3','CDKN3','CENPE','CENPN','DLGAP5','EXO1','GAS6','HMMR','KIF2C','KIF4A','MELK','NCAPD2','NUF2','PBK','RACGAP1','RFC4','TK1','UBE2C','ZWINT']
AR_SIG = ['AR','KLK3','KLK2','TMPRSS2','FKBP5','STEAP2','ACPP','CAMKK2']
EMT = ['CDH1','VIM','CDH2','SNAI1','SNAI2','ZEB1','FN1','CD44','ITGA6']
PROLIF = ['MKI67','TOP2A','PCNA','MCM2','MCM5','MCM7','AURKA','BIRC5','CCNB1']
DNA_REPAIR = ['BRCA1','BRCA2','RAD51','ATM','CHEK2','XRCC2','PARP1','PALB2','RAD54L','GEN1']
CELL_CYCLE = ['CCND1','CCNE1','CDK2','CDK4','CDK6','RB1','E2F1','TP53','CDKN1A','CDKN2A','CDKN1B']
BCR_LIT = ['Prolaris', 'Decipher', 'Cell_Cycle', 'DNA_Repair', 'Proliferation']
SET_MAP = {'Prolaris': PROLARIS, 'Decipher': DECIPHER, 'Cell_Cycle': CELL_CYCLE,
           'DNA_Repair': DNA_REPAIR, 'Proliferation': PROLIF}

def pscores(df):
    s = pd.DataFrame(index=df.index)
    for name, genes in SET_MAP.items():
        p = [g for g in genes if g in df.columns]
        s[name] = df[p].mean(axis=1)
        print(f'  {name}: {len(p)}/{len(genes)} genes present')
    return s

print('TCGA:'); ps_t = pscores(X_train)
print('GSE54460:'); ps_5 = pscores(X54460)

# train on TCGA (C=0.1 as tuned internally on GSE70769 analysis; same config)
lr = LogisticRegression(C=0.1, max_iter=2000, class_weight='balanced')
lr.fit(ps_t, y_train)
p5 = lr.predict_proba(ps_5)[:, 1]
auc5 = roc_auc_score(y54460, p5)

rng = np.random.RandomState(42)
yv = y54460.values
boots = [roc_auc_score(yv[idx], p5[idx])
         for idx in (rng.choice(len(yv), len(yv), replace=True) for _ in range(3000))
         if len(np.unique(yv[idx])) == 2]
lo, hi = np.percentile(boots, [2.5, 97.5])

print(f'\nGSE54460 (RNA-Seq FFPE, n={len(yv)}, BCR+={int(yv.sum())}):')
print(f'  EXTERNAL AUC = {auc5:.3f}  CI [{lo:.3f}, {hi:.3f}]')
print(f'  CI>0.5: {"YES" if lo > 0.5 else "NO"}')

pd.DataFrame([dict(cohort='GSE54460', auc=auc5, ci_lo=lo, ci_hi=hi)]).to_csv(
    'outputs/tables/pathway_gse54460.csv', index=False)
print('Saved -> outputs/tables/pathway_gse54460.csv')
