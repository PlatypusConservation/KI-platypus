#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reproduce_from_deposit.py
=========================
Recomputes the manuscript's headline statistics from the deposited genotype
matrix alone and checks them against the deposited result tables.

    python code/reproduce_from_deposit.py

Run it from the top of the deposit (the folder that contains data/, results/
and code/), or give that folder as the single argument.  It needs nothing but
the files in this deposit: no DArT report, no network, no R.

Inputs
    data/genotypes_filtered_4002x222.csv
    data/samples_analysed_222.csv
Checked against
    data/diversity_by_population.csv
    data/diversity_by_individual.csv
    results/fst_pairwise.csv
    results/pca_eigen.json
    results/pca_scores.csv
    results/sfs.json

Exit status 0 means every check passed.  Each check prints the largest
discrepancy it found, so a failure tells you which statistic moved and by how
much.  The estimators here are written out from the published formulae (Nei
1987; Weir & Cockerham 1984; El Mousadik & Petit 1996; Patterson et al. 2006);
they are the same estimators ki_coprocessed_joint.py uses, applied to the
deposited matrix rather than to the raw DArT report.

Dependencies: python3, numpy, pandas.  Runs in about 20 s.
"""
import itertools
import json
import math
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else '.')
if not os.path.isdir(os.path.join(ROOT, 'data')):
    ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
DATA = os.path.join(ROOT, 'data')
RES = os.path.join(ROOT, 'results')

FAILURES = []
NL = 4002
POPS_ORDER = ['KI', 'EUCUMBENE_ABOVE', 'EUCUMBENE_BELOW', 'MITTA_ABOVE', 'MITTA_BELOW',
              'OVENS', 'SEVERN_ABOVE', 'SEVERN_BELOW', 'SNOWY', 'TENTERFIELD', 'THREDBO']


def report(name, worst, tol, detail=''):
    ok = np.isfinite(worst) and worst <= tol
    print(f'  [{"PASS" if ok else "FAIL"}] {name}: max discrepancy {worst:.3g} (tolerance {tol:g})'
          f'{" " + detail if detail else ""}')
    if not ok:
        FAILURES.append(name)


def hard(name, ok, detail):
    print(f'  [{"PASS" if ok else "FAIL"}] {name}: {detail}')
    if not ok:
        FAILURES.append(name)


# --------------------------------------------------------------- load
print('loading the deposited matrix')
GT = pd.read_csv(os.path.join(DATA, 'genotypes_filtered_4002x222.csv'))
SAMP = pd.read_csv(os.path.join(DATA, 'samples_analysed_222.csv'))

sample_cols = [c for c in GT.columns if c != 'AlleleID']
G = GT[sample_cols].to_numpy(dtype=np.int8)          # loci x individuals

hard('matrix shape 4,002 loci x 222 individuals', G.shape == (NL, 222), f'{G.shape}')
hard('no missing genotypes', int(np.isnan(GT[sample_cols].to_numpy(float)).sum()) == 0
     and bool(((G >= 0) & (G <= 2)).all()), 'every cell is 0, 1 or 2')
hard('sample columns match samples_analysed_222.csv',
     sample_cols == SAMP['sample'].astype(str).tolist(), 'same order, same names')
ki = sorted([c for c in sample_cols if c.startswith('KI')], key=lambda s: int(s[2:]))
hard('KI columns are exactly KI1-KI8', ki == [f'KI{i}' for i in range(1, 9)], ', '.join(ki))
hard('locus keys unique', GT['AlleleID'].is_unique, f'{GT["AlleleID"].nunique()} distinct AlleleIDs')

popv = SAMP['population'].astype(str).to_numpy()
POP_COLS = {p: np.where(popv == p)[0] for p in POPS_ORDER}
hard('population sizes', all(len(POP_COLS[p]) > 0 for p in POPS_ORDER),
     ', '.join(f'{p}={len(POP_COLS[p])}' for p in POPS_ORDER))

_LG = {}


def lgtab(N):
    if N not in _LG:
        _LG[N] = np.array([math.lgamma(k + 1) for k in range(N + 1)])
    return _LG[N]


def locus_stats(dose):
    """Per-locus statistics for one population.  dose: loci x n, no missing."""
    n = dose.shape[1]
    alt = dose.astype(np.float64).sum(1)
    p = alt / (2.0 * n)
    q = 1.0 - p
    Ho = (dose == 1).sum(1) / float(n)
    He_hw = 2.0 * p * q
    He_unb = (2.0 * n / (2.0 * n - 1.0)) * He_hw
    # Weir & Cockerham (1984) eqn 2 with r = 1
    b = (n / (n - 1.0)) * (p * q - ((2.0 * n - 1.0) / (4.0 * n)) * Ho)
    c = Ho / 2.0
    return dict(n=n, p=p, q=q, Ho=Ho, He_hw=He_hw, He_unb=He_unb, b=b, c=c,
                alt_count=np.rint(alt).astype(np.int64))


def allelic_richness(alt_counts, twoN, g=16):
    """El Mousadik & Petit (1996) rarefied allelic richness, biallelic."""
    if twoN < g:
        return float('nan')
    lg = lgtab(twoN)
    lCNg = lg[twoN] - lg[g] - lg[twoN - g]
    nref = twoN - alt_counts
    r_alt = np.where(nref >= g, np.exp((lg[nref] - lg[g] - lg[np.maximum(nref - g, 0)]) - lCNg), 0.0)
    r_ref = np.where(alt_counts >= g,
                     np.exp((lg[alt_counts] - lg[g] - lg[np.maximum(alt_counts - g, 0)]) - lCNg), 0.0)
    return float(((1.0 - r_alt) + (1.0 - r_ref)).mean())


ST = {p: locus_stats(G[:, POP_COLS[p]]) for p in POPS_ORDER}
ST_JOINT = locus_stats(G)

# --------------------------------------------------------------- per population
print('\nper-population diversity (data/diversity_by_population.csv)')
POP = pd.read_csv(os.path.join(DATA, 'diversity_by_population.csv')).set_index('population')
worst = {k: 0.0 for k in ['Ho', 'He_unbiased', 'Fis_WC', 'Ar_g16', 'mean_folded_MAF']}
npoly_bad = []
for p in POPS_ORDER:
    s = ST[p]
    row = POP.loc[p]
    got = {
        'Ho': float(s['Ho'].mean()),
        'He_unbiased': float(s['He_unb'].mean()),
        'Fis_WC': float(s['b'].sum() / (s['b'] + s['c']).sum()),
        'Ar_g16': allelic_richness(s['alt_count'], 2 * s['n'], 16),
        'mean_folded_MAF': float(np.minimum(s['p'], s['q'])[(s['p'] > 0) & (s['p'] < 1)].mean()),
    }
    for k, v in got.items():
        ref = float(row[k])
        if np.isnan(ref) and np.isnan(v):
            continue
        worst[k] = max(worst[k], abs(v - ref))
    npoly = int(((s['p'] > 0) & (s['p'] < 1)).sum())
    if npoly != int(row['n_polymorphic']):
        npoly_bad.append((p, npoly, int(row['n_polymorphic'])))
for k, v in worst.items():
    report(f'{k} over 11 populations', v, 1e-9)
hard('polymorphic-locus counts over 11 populations', not npoly_bad,
     'all 11 reproduced exactly' if not npoly_bad else str(npoly_bad))

# --------------------------------------------------------------- per individual
print('\nper-individual heterozygosity (data/diversity_by_individual.csv)')
IND = pd.read_csv(os.path.join(DATA, 'diversity_by_individual.csv')).set_index('sample')
ho = (G == 1).mean(axis=0)
report('H_O over 222 individuals', float(np.max(np.abs(ho - IND.loc[sample_cols, 'Ho_ind'].values))),
       1e-9)

# --------------------------------------------------------------- Fst
print('\npairwise differentiation (results/fst_pairwise.csv)')


def wc_components(s1, s2):
    n1, n2, r = float(s1['n']), float(s2['n']), 2.0
    nbar = (n1 + n2) / r
    nc = (r * nbar - (n1 * n1 + n2 * n2) / (r * nbar)) / (r - 1.0)
    pbar = (n1 * s1['p'] + n2 * s2['p']) / (n1 + n2)
    s2_ = (n1 * (s1['p'] - pbar) ** 2 + n2 * (s2['p'] - pbar) ** 2) / ((r - 1.0) * nbar)
    hbar = (n1 * s1['Ho'] + n2 * s2['Ho']) / (n1 + n2)
    core = pbar * (1.0 - pbar) - ((r - 1.0) / r) * s2_
    a = (nbar / nc) * (s2_ - (1.0 / (nbar - 1.0)) * (core - hbar / 4.0))
    b = (nbar / (nbar - 1.0)) * (core - ((2.0 * nbar - 1.0) / (4.0 * nbar)) * hbar)
    return a, b, hbar / 2.0


FST = pd.read_csv(os.path.join(RES, 'fst_pairwise.csv')).set_index(['pop_i', 'pop_j'])
wf = 0.0
for pi, pj in itertools.combinations(POPS_ORDER, 2):
    a, b, c = wc_components(ST[pi], ST[pj])
    got = float(a.sum() / (a + b + c).sum())
    key = (pi, pj) if (pi, pj) in FST.index else (pj, pi)
    wf = max(wf, abs(got - float(FST.loc[key, 'Fst_WC'])))
report('F_ST (Weir & Cockerham) over all 55 pairs', wf, 1e-9)

FMAT = pd.read_csv(os.path.join(RES, 'fst_matrix.csv')).set_index('population')
wm = 0.0
for pi, pj in itertools.combinations(POPS_ORDER, 2):
    key = (pi, pj) if (pi, pj) in FST.index else (pj, pi)
    wm = max(wm, abs(float(FMAT.loc[pi, pj]) - float(FST.loc[key, 'Fst_WC'])),
             abs(float(FMAT.loc[pj, pi]) - float(FST.loc[key, 'Fst_WC'])))
report('fst_matrix.csv agrees with fst_pairwise.csv', wm, 1e-12)

# --------------------------------------------------------------- PCA
print('\nprincipal component analysis (results/pca_eigen.json, results/pca_scores.csv)')
pj_ = ST_JOINT['p']
qj_ = 1.0 - pj_
scale = np.where(pj_ * qj_ > 0, np.sqrt(2.0 * pj_ * qj_), 1.0)
X = (G.astype(np.float64).T - 2.0 * pj_) / scale
U, S, _ = np.linalg.svd(X, full_matrices=False)
scores = U * S
varpct = 100.0 * (S ** 2) / float((S ** 2).sum())
EIG = json.loads(open(os.path.join(RES, 'pca_eigen.json')).read())
report('variance explained, PC1-PC10',
       float(np.max(np.abs(varpct[:10] - np.asarray(EIG['variance_pct'][:10])))), 1e-8)

PCS = pd.read_csv(os.path.join(RES, 'pca_scores.csv'))
hard('pca_scores.csv sample order', PCS['sample'].astype(str).tolist() == sample_cols,
     'matches the genotype matrix column order')
wp = 0.0
for k in range(6):
    ref = PCS[f'PC{k + 1}'].to_numpy(float)
    got = scores[:, k]
    wp = max(wp, float(min(np.max(np.abs(got - ref)), np.max(np.abs(got + ref)))))
report('PC1-PC6 scores (sign-invariant)', wp, 1e-6)

is_ki = np.isin(sample_cols, ki)
pc2 = scores[:, 1] * (1.0 if abs(scores[is_ki, 1].mean() - PCS.loc[is_ki, 'PC2'].mean()) <
                      abs(-scores[is_ki, 1].mean() - PCS.loc[is_ki, 'PC2'].mean()) else -1.0)
hard('KI fully separated from all 214 mainland individuals on PC2',
     pc2[is_ki].min() > pc2[~is_ki].max() or pc2[is_ki].max() < pc2[~is_ki].min(),
     f'KI range [{pc2[is_ki].min():+.2f}, {pc2[is_ki].max():+.2f}] against mainland '
     f'[{pc2[~is_ki].min():+.2f}, {pc2[~is_ki].max():+.2f}]')

# --------------------------------------------------------------- SFS
print('\nKI folded site-frequency spectrum (results/sfs.json)')
SFS = json.loads(open(os.path.join(RES, 'sfs.json')).read())
sKI = ST['KI']
poly = (sKI['p'] > 0) & (sKI['p'] < 1)
mac = np.minimum(sKI['alt_count'], 16 - sKI['alt_count'])[poly]
counts = {str(k): int((mac == k).sum()) for k in range(1, 9)}
hard('KI polymorphic loci', int(poly.sum()) == int(SFS['KI']['n_polymorphic']),
     f'{int(poly.sum())} (deposited {SFS["KI"]["n_polymorphic"]})')
hard('KI minor-allele-count spectrum', counts == {k: int(v) for k, v in SFS['KI']['mac_counts'].items()},
     ', '.join(f'{k}:{v}' for k, v in counts.items()))
report('KI mean folded MAF',
       abs(float(np.minimum(sKI['p'], sKI['q'])[poly].mean()) - float(SFS['KI']['mean_folded_maf'])),
       1e-12)

OV = ST['OVENS']
povp = (OV['p'] > 0) & (OV['p'] < 1)
report('Ovens full-sample (n = 19) mean folded MAF',
       abs(float(np.minimum(OV['p'], OV['q'])[povp].mean()) -
           float(SFS['OVENS_full']['mean_folded_maf'])), 1e-12,
       f'({int(povp.sum())} polymorphic loci)')

# --------------------------------------------------------------- abstract numbers
print('\nabstract and Table 1 headline values')
he = {p: float(ST[p]['He_unb'].mean()) for p in POPS_ORDER}
ml_mean = float(np.mean([he[p] for p in POPS_ORDER if p != 'KI']))
hard('KI H_E = 0.0954', abs(he['KI'] - 0.0954) < 5e-5, f'{he["KI"]:.4f}')
hard('KI H_E is 86.5% of the mainland mean',
     abs(100.0 * he['KI'] / ml_mean - 86.5) < 0.05,
     f'{100.0 * he["KI"] / ml_mean:.1f}% of {ml_mean:.4f}')
hard('KI is the lowest of the eleven populations for H_E',
     min(he, key=he.get) == 'KI', f'lowest = {min(he, key=he.get)}')
ar = {p: allelic_richness(ST[p]['alt_count'], 2 * ST[p]['n'], 16) for p in POPS_ORDER}
ar_ok = {p: v for p, v in ar.items() if np.isfinite(v)}
hard('KI has the lowest A_r(g = 16) of the nine standardisable populations',
     min(ar_ok, key=ar_ok.get) == 'KI' and len(ar_ok) == 9,
     f'{len(ar_ok)} populations standardisable, lowest = {min(ar_ok, key=ar_ok.get)} '
     f'at {ar_ok["KI"]:.4f}')
kifst = {}
for p in POPS_ORDER:
    if p == 'KI':
        continue
    key = ('KI', p) if ('KI', p) in FST.index else (p, 'KI')
    kifst[p] = float(FST.loc[key, 'Fst_WC'])
hard('KI F_ST range 0.147-0.264, mean 0.191',
     abs(min(kifst.values()) - 0.1468) < 5e-4 and abs(max(kifst.values()) - 0.2637) < 5e-4
     and abs(np.mean(list(kifst.values())) - 0.1911) < 5e-4,
     f'{min(kifst.values()):.4f} ({min(kifst, key=kifst.get)}) to {max(kifst.values()):.4f} '
     f'({max(kifst, key=kifst.get)}), mean {np.mean(list(kifst.values())):.4f}')
fis_ki = float(ST['KI']['b'].sum() / (ST['KI']['b'] + ST['KI']['c']).sum())
hard('KI F_IS = -0.0385', abs(fis_ki + 0.0385) < 5e-5, f'{fis_ki:+.4f}')

# --------------------------------------------------------------- verdict
print()
if FAILURES:
    print(f'{len(FAILURES)} CHECK(S) FAILED: ' + '; '.join(FAILURES))
    sys.exit(1)
print('all checks passed: the deposited result tables are reproduced from the '
      'deposited genotype matrix')
sys.exit(0)
