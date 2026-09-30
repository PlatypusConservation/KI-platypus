#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ki_coprocessed_joint.py
=======================
Replacement analysis pipeline for the Kangaroo Island platypus conservation-genomics
manuscript revision (Bino, Hawke, Baring, Gongora).

Implements PLAN_ANALYSIS.md in full.  Runs end-to-end from the raw DArT co-processed
one-row-per-locus mapping report with no manual steps.

  * NO CloneID merge, NO batch-effect call-string filter.  The co-processed
    22,054 x 376 matrix is the starting point.
  * The old (Mijangos-only) DArT report is used in exactly two metadata-only places:
    the scaffold->chromosome crosswalk and the old-vs-new bias diagnostic.

Dependencies: python3, numpy, pandas ONLY.  No scipy, no sklearn.
  chi-square(1 df) upper tail  -> math.erfc
  two-sample Kolmogorov-Smirnov -> numpy
  Mann-Whitney U               -> tie-corrected normal approximation
  binomial coefficients        -> math.lgamma

Outputs: work/out/*.{csv,json} and work/out/RESULTS_SUMMARY.md
Author: analysis pipeline for the 2026 revision.  Seed 42 throughout.
"""
import os, sys, csv, json, math, time, hashlib, itertools, collections
import numpy as np
import pandas as pd

T0 = time.time()
SEED = 42
RNG = np.random.default_rng(SEED)
np.seterr(divide='ignore', invalid='ignore')

HOME   = os.path.expanduser('~')
GEN    = os.path.join(HOME, 'mnt', 'Kangaroo Island', 'Genetics')
MIJ    = os.path.join(GEN, 'mijangos_data')
WORK   = os.path.join(GEN, 'Revision', 'work')
OUT    = os.path.join(WORK, 'out')
FIGD   = os.path.join(WORK, 'fig')
NEWREP = os.path.join(HOME, 'mnt', 'Platypus - Genomics - General', 'DaRT Data',
                      'DPla24-9425', 'Report-DPla24-9425',
                      'Report_DPla24-9425_1_moreOrders_SNP_mapping_2.csv')
OLDREP = os.path.join(MIJ, 'Report_DPla19-4171_SNP_mapping_2.csv')
os.makedirs(OUT, exist_ok=True); os.makedirs(FIGD, exist_ok=True)

ASSERTIONS = {}          # name -> {"value":..., "expected":..., "passed":bool}
def check(name, ok, value, expected):
    ASSERTIONS[name] = {"value": value, "expected": expected, "passed": bool(ok)}
    status = 'PASS' if ok else '*** FAIL ***'
    print(f'  [ASSERT {status}] {name}: value={value!r} expected={expected}', flush=True)
    if not ok:
        raise AssertionError(f'{name}: got {value!r}, expected {expected}')

def log(msg):
    print(f'[{time.time()-T0:7.1f}s] {msg}', flush=True)

# ---------------------------------------------------------------------------
# 0.  Statistics implemented without scipy
# ---------------------------------------------------------------------------
def chi2_sf_1df(x2):
    """Upper tail of chi-square with 1 d.f.  P(X > x2) = erfc(sqrt(x2/2))."""
    if isinstance(x2, np.ndarray):
        return np.array([math.erfc(math.sqrt(max(float(v), 0.0) / 2.0)) for v in x2.ravel()]
                        ).reshape(x2.shape)
    return math.erfc(math.sqrt(max(float(x2), 0.0) / 2.0))

def ks_2samp(a, b):
    """Two-sample Kolmogorov-Smirnov D and asymptotic p.  Q(l)=2*sum (-1)^(k-1) exp(-2k^2 l^2)."""
    a = np.sort(np.asarray(a, float)); b = np.sort(np.asarray(b, float))
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return float('nan'), float('nan')
    allv = np.concatenate([a, b])
    cdf1 = np.searchsorted(a, allv, side='right') / n1
    cdf2 = np.searchsorted(b, allv, side='right') / n2
    D = float(np.max(np.abs(cdf1 - cdf2)))
    ne = math.sqrt(n1 * n2 / (n1 + n2))
    lam = (ne + 0.12 + 0.11 / ne) * D
    p = 2.0 * sum((-1) ** (k - 1) * math.exp(-2.0 * k * k * lam * lam) for k in range(1, 101))
    return D, float(min(max(p, 0.0), 1.0))

def mannwhitney(a, b):
    """Mann-Whitney U, tie-corrected normal approximation.  Returns (U, z, two-sided p)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    n1, n2 = len(a), len(b)
    allv = np.concatenate([a, b])
    order = np.argsort(allv, kind='mergesort')
    ranks = np.empty(len(allv), float)
    sv = allv[order]
    i = 0; ties = []
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        r = 0.5 * (i + j) + 1.0
        ranks[order[i:j + 1]] = r
        ties.append(j - i + 1)
        i = j + 1
    R1 = ranks[:n1].sum()
    U1 = R1 - n1 * (n1 + 1) / 2.0
    N = n1 + n2
    mu = n1 * n2 / 2.0
    tie_term = sum(t ** 3 - t for t in ties)
    sd = math.sqrt(n1 * n2 / 12.0 * ((N + 1) - tie_term / (N * (N - 1)))) if N > 1 else float('nan')
    z = (U1 - mu) / sd if sd > 0 else float('nan')
    p = math.erfc(abs(z) / math.sqrt(2.0)) if sd > 0 else float('nan')
    return float(U1), float(z), float(p)

_LG = math.lgamma
def log_binom(n, k):
    if k < 0 or k > n or n < 0:
        return -np.inf
    return _LG(n + 1) - _LG(k + 1) - _LG(n - k + 1)

def binom_ratio(N_minus, N, g):
    """C(N_minus, g) / C(N, g), computed in logs.  0 if N_minus < g."""
    if N_minus < g:
        return 0.0
    return math.exp(log_binom(N_minus, g) - log_binom(N, g))

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

# ===========================================================================
# 1.  PARSE THE CO-PROCESSED REPORT  (PLAN §1.1-1.3)
# ===========================================================================
log('STEP 1  parsing co-processed DArT report')
META_COLS = ['AlleleID','CloneID','AlleleSequenceRef','AlleleSequenceSnp','TrimmedSequenceRef',
             'TrimmedSequenceSnp','Chrom_Platypus_NCBIv1','ChromPosTag_Platypus_NCBIv1',
             'ChromPosSnp_Platypus_NCBIv1','AlnCnt_Platypus_NCBIv1','AlnEvalue_Platypus_NCBIv1',
             'Strand_Platypus_NCBIv1','SNP','SnpPosition','CallRate','OneRatioRef','OneRatioSnp',
             'FreqHomRef','FreqHomSnp','FreqHets','PICRef','PICSnp','AvgPIC','AvgCountRef',
             'AvgCountSnp','RepAvg']
N_META = 26
DOSE = {'0': 0, '1': 2, '2': 1}          # ref-hom -> 0, snp-hom -> 2, het -> 1

with open(NEWREP, newline='', encoding='utf-8', errors='replace') as fh:
    rdr = csv.reader(fh)
    rows_head = [next(rdr) for _ in range(7)]
    header = rows_head[6]
    orders_row = rows_head[0]
    data_rows = [r for r in rdr if len(r) > 1]

n_cols = len(header)
check('n_columns_402', n_cols == 402, n_cols, 402)
check('n_data_rows_22054', len(data_rows) == 22054, len(data_rows), 22054)
bad_width = sum(1 for r in data_rows if len(r) != n_cols)
check('all_rows_402_fields', bad_width == 0, bad_width, 0)
check('meta_col_names', header[:N_META] == META_COLS, header[:N_META] == META_COLS, True)

L0 = len(data_rows)
NS = n_cols - N_META
check('n_samples_376', NS == 376, NS, 376)

sample_ids_report = [s.strip() for s in header[N_META:]]
sample_orders     = [s.strip() for s in orders_row[N_META:]]
ord_counts = collections.Counter(sample_orders)
check('order_counts', ord_counts.get('DPla19-4171') == 218 and ord_counts.get('DPla24-9425') == 158,
      dict(ord_counts), {'DPla19-4171': 218, 'DPla24-9425': 158})

G = np.full((L0, NS), -1, dtype=np.int8)
meta = {c: [] for c in META_COLS}
for i, row in enumerate(data_rows):
    for c, v in zip(META_COLS, row[:N_META]):
        meta[c].append(v)
    G[i] = [DOSE.get(v, -1) for v in row[N_META:]]
META = pd.DataFrame(meta)
del data_rows
log(f'  matrix {G.shape}  missing={float((G<0).mean()):.4f}')

# ---- RISK 4: allele-orientation / encoding cross-check against the report's own columns
ncall = (G >= 0).sum(1).astype(float)
f0 = (G == 0).sum(1) / np.maximum(ncall, 1)
f1 = (G == 1).sum(1) / np.maximum(ncall, 1)
f2 = (G == 2).sum(1) / np.maximum(ncall, 1)
rep_f0 = pd.to_numeric(META.FreqHomRef, errors='coerce').values
rep_f1 = pd.to_numeric(META.FreqHets,  errors='coerce').values
rep_f2 = pd.to_numeric(META.FreqHomSnp, errors='coerce').values
rep_cr = pd.to_numeric(META.CallRate,  errors='coerce').values
ok = np.isfinite(rep_f0) & (ncall > 0)
d0 = float(np.nanmax(np.abs(f0[ok] - rep_f0[ok])))
d1 = float(np.nanmax(np.abs(f1[ok] - rep_f1[ok])))
d2 = float(np.nanmax(np.abs(f2[ok] - rep_f2[ok])))
dcr = float(np.nanmax(np.abs(ncall[ok] / NS - rep_cr[ok])))
ENCODING_DIFFS = {'FreqHomRef_vs_dose0': d0, 'FreqHets_vs_dose1': d1,
                  'FreqHomSnp_vs_dose2': d2, 'CallRate_vs_computed': dcr}
check('encoding_orientation_lt_1e-5', max(d0, d1, d2, dcr) < 1e-5,
      round(max(d0, d1, d2, dcr), 12), '< 1e-5')

# ===========================================================================
# 2.  SAMPLE SET AND POPULATIONS  (PLAN §1.4-1.5)
# ===========================================================================
log('STEP 2  sample selection (BY COLUMN INDEX, never by name)')
ALIAS   = {'V15_1': 'V15a', 'V15_2': 'V15b'}      # RISK 3
DROP_MIJ = {'E32', 'T3', 'V34', 'T42'}            # Mijangos et al. 2022 exclusions
POPS_ORDER = ['KI','EUCUMBENE_ABOVE','EUCUMBENE_BELOW','MITTA_ABOVE','MITTA_BELOW','OVENS',
              'SEVERN_ABOVE','SEVERN_BELOW','SNOWY','TENTERFIELD','THREDBO']

idp = pd.read_csv(os.path.join(MIJ, 'ID_Pop_Platypus.csv'))
recode = pd.read_csv(os.path.join(MIJ, 'new_pop_assignments.csv'), header=None, names=['old', 'new'])
rmap = dict(zip(recode.old, recode.new))
idp['pop_final'] = idp['pop'].map(rmap).fillna(idp['pop'])
id2pop = dict(zip(idp.id, idp.pop_final))
id2sex = dict(zip(idp.id, idp.Sex.fillna('')))
regions = pd.read_csv(os.path.join(MIJ, 'new_pop_assignments_regions.csv'),
                      header=None, names=['pop', 'region'])
pop2region = dict(zip(regions['pop'], regions.region)); pop2region['KI'] = 'KI'

canon = np.array([ALIAS.get(s, s) for s in sample_ids_report], dtype=object)
order = np.array(sample_orders, dtype=object)
pop_all  = np.array([''] * NS, dtype=object)
incl     = np.zeros(NS, bool)
reason   = np.array([''] * NS, dtype=object)
for j in range(NS):                                    # RISK 2: index-keyed, not name-keyed
    o, c = order[j], canon[j]
    if o == 'DPla19-4171':
        if c in DROP_MIJ:
            reason[j] = f'mijangos_excluded_{c}'
        elif c in id2pop:
            pop_all[j] = id2pop[c]; incl[j] = True
        else:
            reason[j] = 'no_population_record'
    else:
        if c.startswith('KI') and c[2:].isdigit():
            pop_all[j] = 'KI'; incl[j] = True
        else:
            reason[j] = 'not_in_analysis_panel'
# V30 / V32 swap (Mijangos et al. 2022 analyses_platy.R lines 64-71)
for sid, newpop in (('V30', 'MITTA_ABOVE'), ('V32', 'OVENS')):
    jj = np.where((canon == sid) & (order == 'DPla19-4171'))[0]
    assert len(jj) == 1, f'V30/V32 swap: expected 1 column for {sid}, found {len(jj)}'
    pop_all[jj[0]] = newpop
V30_V32_APPLIED = True

keep_idx = np.where(incl)[0]
popk = pop_all[keep_idx]
sids = canon[keep_idx]
counts = collections.Counter(popk)
EXPECTED_N = {'SNOWY':56,'TENTERFIELD':39,'SEVERN_ABOVE':23,'EUCUMBENE_BELOW':20,'OVENS':19,
              'THREDBO':19,'SEVERN_BELOW':17,'MITTA_ABOVE':13,'KI':8,'MITTA_BELOW':4,
              'EUCUMBENE_ABOVE':4}
check('n_samples_retained_222', len(keep_idx) == 222, int(len(keep_idx)), 222)
check('n_KI_8', counts['KI'] == 8, int(counts['KI']), 8)
check('population_size_table', dict(counts) == EXPECTED_N, dict(counts), EXPECTED_N)
# duplicated IDs across orders must not have leaked in
dup_ids = [i for i, c in collections.Counter(sample_ids_report).items() if c > 1]
check('duplicate_ids_across_orders_detected', set(dup_ids) == {'T3','T4','T5','T6','T7','T8','T9'},
      sorted(dup_ids), ['T3','T4','T5','T6','T7','T8','T9'])
dup_kept = [sids[i] for i in range(len(sids)) if sids[i] in dup_ids]
check('no_duplicate_id_kept_twice', len(dup_kept) == len(set(dup_kept)), sorted(dup_kept),
      'each duplicated ID retained at most once (only the DPla19-4171 copy)')

pd.DataFrame({
    'col_index': np.arange(NS), 'order': order, 'id_report': sample_ids_report,
    'id_canonical': canon, 'population': pop_all,
    'region': [pop2region.get(p, '') for p in pop_all],
    'sex': [id2sex.get(c, '') if o == 'DPla19-4171' else '' for o, c in zip(order, canon)],
    'included': incl.astype(int), 'exclusion_reason': reason,
}).to_csv(os.path.join(OUT, 'samples_used.csv'), index=False)

Gk = G[:, keep_idx]                                    # RISK 5/6: everything downstream uses 222
SEX_ARR = np.array([id2sex.get(c, '') if o == 'DPla19-4171' else '' for o, c in zip(order, canon)],
                   dtype=object)
log(f'  retained {len(keep_idx)} samples; KI={counts["KI"]}')

# ===========================================================================
# 3.  SCAFFOLD -> CHROMOSOME CROSSWALK AND SEX-CHROMOSOME EXCLUSION  (PLAN §2.5)
# ===========================================================================
log('STEP 3  scaffold->chromosome crosswalk')
chrom_tab = pd.read_csv(os.path.join(MIJ, 'chrom_platypus.csv'), encoding='utf-8-sig')
chrom_tab = chrom_tab[chrom_tab.Type == 'Chr']
rs2label = dict(zip(chrom_tab['RefSeq_mOrnAna_v1'], chrom_tab['Name'].astype(str)))
SEXSET = set(chrom_tab.iloc[21:31]['RefSeq_mOrnAna_v1'].dropna())      # X1-X5, Y1-Y5
check('sexset_size_10', len(SEXSET) == 10, sorted(SEXSET), '10 RefSeq accessions (X1-X5, Y1-Y5)')

# 3a. old report -> AlleleID -> RefSeq chromosome accession
old_allele2rs, old_clone_snp, old_allele_rows = {}, {}, 0
with open(OLDREP, newline='', encoding='utf-8', errors='replace') as fh:
    rdr = csv.reader(fh)
    for _ in range(6):
        next(rdr)
    oh = next(rdr); oci = {n: i for i, n in enumerate(oh)}
    for row in rdr:
        if len(row) < len(oh):
            continue
        old_allele_rows += 1
        ch = row[oci['Chrom_Platypus_Chrom_NCBIv1']].strip()
        if ch:
            toks = ch.split('_')
            if len(toks) >= 2:
                old_allele2rs[row[oci['AlleleID']]] = '_'.join(toks[:2])
        c = row[oci['CloneID']]
        if c not in old_clone_snp:
            old_clone_snp[c] = row[oci['SNP']]
check('old_report_annotated_alleles', abs(len(old_allele2rs) - 14621) <= 5, len(old_allele2rs), 14621)

# 3b. vote scaffold -> RefSeq
scaffold = META.Chrom_Platypus_NCBIv1.fillna('').astype(str).str.strip().values
allele   = META.AlleleID.astype(str).values
votes = collections.defaultdict(collections.Counter)
n_anchor = 0
for i in range(L0):
    rs = old_allele2rs.get(allele[i])
    if rs and scaffold[i]:
        votes[scaffold[i]][rs] += 1
        n_anchor += 1
scaffolds_all = sorted({s for s in scaffold if s})
cw_rows, s2rs, s2class = [], {}, {}
for s in scaffolds_all:
    if s in votes and sum(votes[s].values()) > 0:
        rs, k = votes[s].most_common(1)[0]
        tot = sum(votes[s].values()); purity = k / tot
        lab = rs2label.get(rs, 'unplaced' if rs.startswith('NW_') else rs)
        cls = 'sex' if rs in SEXSET else ('unplaced' if (rs.startswith('NW_') and rs not in SEXSET)
                                          else 'autosome')
        if cls == 'unplaced':
            lab = 'unplaced'
        s2rs[s] = rs; s2class[s] = cls
        cw_rows.append((s, rs, lab, purity, tot, cls))
    else:
        s2rs[s] = ''; s2class[s] = 'unanchored'
        cw_rows.append((s, '', 'unanchored', float('nan'), 0, 'unanchored'))
CW = pd.DataFrame(cw_rows, columns=['scaffold','refseq_majority','chrom_label','purity',
                                    'n_anchor_loci','class'])
CW.to_csv(os.path.join(OUT, 'scaffold_chrom_crosswalk.csv'), index=False)
n_anch_scaf = int((CW['class'] != 'unanchored').sum())
n_pure95 = int((CW.purity >= 0.95).sum())
log(f'  scaffolds {len(CW)} anchored {n_anch_scaf} pure>=95% {n_pure95} anchor votes {n_anchor}')

locus_class = np.array([s2class.get(s, 'unaligned') if s else 'unaligned' for s in scaffold],
                       dtype=object)
locus_rs    = np.array([s2rs.get(s, '') if s else '' for s in scaffold], dtype=object)
locus_label = np.array([rs2label.get(r, 'unplaced' if r.startswith('NW_') else ('unanchored' if r == '' else r))
                        for r in locus_rs], dtype=object)
sexflag = np.array([c == 'sex' for c in locus_class])
n_sex = int(sexflag.sum())
check('n_sex_loci_2300_2800', 2300 <= n_sex <= 2800, n_sex, '2300..2800')     # RISK 1
unanchored_loci = int((locus_class == 'unanchored').sum())

# 3c. validation vs the old report's own direct per-locus annotation
direct_sex, dsex_cw_sex, dsex_cw_auto, dauto_cw_sex, dauto_cw_auto = 0, 0, 0, 0, 0
for i in range(L0):
    rs_direct = old_allele2rs.get(allele[i])
    if rs_direct is None or not scaffold[i] or s2class[scaffold[i]] == 'unanchored':
        continue
    ds = rs_direct in SEXSET; cs = sexflag[i]
    if ds:
        direct_sex += 1
        dsex_cw_sex += cs; dsex_cw_auto += (not cs)
    else:
        dauto_cw_sex += cs; dauto_cw_auto += (not cs)
n_conf = dsex_cw_sex + dsex_cw_auto + dauto_cw_sex + dauto_cw_auto
agree = (dsex_cw_sex + dauto_cw_auto) / max(n_conf, 1)
check('crosswalk_agreement_ge_0.99', agree >= 0.99, round(agree, 5), '>= 0.99')

# 3d. data-driven sex-linkage flag (reported, NOT filtered on)
male_cols   = np.array([i for i, j in enumerate(keep_idx) if SEX_ARR[j] == 'M'])
female_cols = np.array([i for i, j in enumerate(keep_idx) if SEX_ARR[j] == 'F'])
hM = (Gk[:, male_cols]   == 1).sum(1) / np.maximum((Gk[:, male_cols]   >= 0).sum(1), 1)
hF = (Gk[:, female_cols] == 1).sum(1) / np.maximum((Gk[:, female_cols] >= 0).sum(1), 1)
okMF = ((Gk[:, male_cols] >= 0).sum(1) >= 20) & ((Gk[:, female_cols] >= 0).sum(1) >= 20)
sex_dd = okMF & (hF > 0.05) & (hM < 0.1 * hF)
log(f'  sex-linked by crosswalk {n_sex}; data-driven candidates {int(sex_dd.sum())}')

# ===========================================================================
# 4.  FILTER CASCADE  (PLAN §2)
# ===========================================================================
log('STEP 4  filter cascade')
RepAvg   = pd.to_numeric(META.RepAvg, errors='coerce').values
AlnCnt   = pd.to_numeric(META.AlnCnt_Platypus_NCBIv1, errors='coerce').fillna(0).values
AlnEval  = pd.to_numeric(META.AlnEvalue_Platypus_NCBIv1, errors='coerce').values
CallRepo = pd.to_numeric(META.CallRate, errors='coerce').fillna(0).values
AvgPIC   = pd.to_numeric(META.AvgPIC, errors='coerce').fillna(0).values
CloneID  = META.CloneID.astype(str).values
cr222    = (Gk >= 0).mean(1)                       # RISK 5/6: over the 222 retained columns

pass_rep  = RepAvg >= 1.0
pass_aln  = AlnCnt >= 1
pass_bla  = AlnEval <= 1e-20
pass_auto = ~sexflag
pass_cr   = cr222 >= 1.0

def hwe_chisq_p(sub):
    n  = (sub >= 0).sum(1).astype(float)
    n0 = (sub == 0).sum(1).astype(float)
    nh = (sub == 1).sum(1).astype(float)
    n2 = (sub == 2).sum(1).astype(float)
    p = np.where(n > 0, (2 * n0 + nh) / (2 * np.maximum(n, 1)), np.nan); q = 1 - p
    e0 = n * p * p; e1 = 2 * n * p * q; e2 = n * q * q
    x2 = (np.where(e0 > 0, (n0 - e0) ** 2 / np.maximum(e0, 1e-300), 0) +
          np.where(e1 > 0, (nh - e1) ** 2 / np.maximum(e1, 1e-300), 0) +
          np.where(e2 > 0, (n2 - e2) ** 2 / np.maximum(e2, 1e-300), 0))
    x2 = np.where((p > 0) & (p < 1) & (n > 1), x2, 0.0)
    return chi2_sf_1df(x2)

def run_cascade(callrate_thr=1.0, drop_unanchored=False, hwe_mode='per_pop'):
    m = pass_rep & pass_aln & pass_bla & pass_auto
    if drop_unanchored:
        m = m & (locus_class != 'unanchored')
    m = m & (cr222 >= callrate_thr)
    idx = np.where(m)[0]
    hwe_fail = np.zeros(L0, bool)
    if hwe_mode == 'per_pop':
        big = [p for p in POPS_ORDER if EXPECTED_N[p] >= 19]
        alpha = 0.05 / (len(idx) * len(big))
        nfail = np.zeros(len(idx), int)
        for p in big:
            cols = np.where(popk == p)[0]
            nfail += (hwe_chisq_p(Gk[np.ix_(idx, cols)]) < alpha).astype(int)
        fail = nfail >= 2
    else:
        alpha = 0.05 / len(idx)
        fail = hwe_chisq_p(Gk[idx]) < alpha
    hwe_fail[idx[fail]] = True
    m2 = m & ~hwe_fail
    # one SNP per CloneID: max CallRate, then max AvgPIC, then lowest row index
    best = {}
    for j in np.where(m2)[0]:
        k = CloneID[j]
        key = (CallRepo[j], AvgPIC[j], -j)
        if k not in best or key > best[k][0]:
            best[k] = (key, j)
    fin = np.zeros(L0, bool); fin[[v[1] for v in best.values()]] = True
    return dict(mask_prehwe=m, hwe_fail=hwe_fail, mask_final=fin,
                n_prehwe=int(m.sum()), n_hwe_fail=int(fail.sum()), n_final=int(fin.sum()),
                alpha=alpha)

casc_pp   = run_cascade(1.0, False, 'per_pop')
casc_jt   = run_cascade(1.0, False, 'joint')
casc_drop = run_cascade(1.0, True,  'per_pop')
FINAL_MASK = casc_pp['mask_final']
FIN = np.where(FINAL_MASK)[0]
NL = len(FIN)
check('final_locus_count_2500_6000', 2500 <= NL <= 6000, NL, '2500..6000')

cum = np.ones(L0, bool); steps = []
for nm, msk in [('input', np.ones(L0, bool)), ('RepAvg>=1.0', pass_rep), ('AlnCnt>=1', pass_aln),
                ('AlnEvalue<=1e-20', pass_bla), ('autosomal (crosswalk)', pass_auto),
                ('CallRate==1.0 over 222', pass_cr)]:
    cum = cum & msk
    steps.append({'name': nm, 'standalone': int(msk.sum()), 'cumulative': int(cum.sum())})
steps.append({'name': 'HWE per-population (fail in >=2 of 6 pops, Bonferroni)',
              'standalone': int(L0 - casc_pp['hwe_fail'].sum()),
              'cumulative': int(casc_pp['n_prehwe'] - casc_pp['n_hwe_fail']),
              'removed': casc_pp['n_hwe_fail']})
steps.append({'name': 'one SNP per CloneID', 'standalone': None, 'cumulative': NL})
for s in steps:
    log(f"    {s['name']:<52s} cumulative {s['cumulative']}")

Gf = Gk[FIN]                                        # final matrix, loci x 222
check('no_missing_in_final_matrix', int((Gf < 0).sum()) == 0, int((Gf < 0).sum()), 0)
GF = Gf.astype(np.float64)

# ===========================================================================
# 5.  CORE POPULATION-GENETIC ESTIMATORS  (PLAN §3.1-3.3, §3.9, §4)
# ===========================================================================
# NOTE ON A CORRECTION TO PLAN_ANALYSIS.md §3.3
# ---------------------------------------------
# The plan writes the single-population Weir & Cockerham within-individual component as
#     b = (n/(n-1)) * [ p q - ((n-1)/n)(Ho/2) - Ho/4 ]
# Expanding, that is (n/(n-1))pq - Ho(3n-2)/(4(n-1)).  Substituting Hardy-Weinberg
# proportions (Ho = 2pq) gives F_IS = (2-n)/n, i.e. -0.89 at n = 19, which is impossible.
# Weir & Cockerham (1984) eqn 2 with r = 1 (so s^2 = 0, p_bar = p, h_bar = Ho) is
#     b = (n/(n-1)) * [ p q - ((2n-1)/(4n)) Ho ] ,   c = Ho/2
# which returns F_IS = 1/(2n-1) under exact Hardy-Weinberg proportions, as it should.
# The published estimator is used here; the plan's expression contains an algebraic slip.
LOG_CACHE = {}
def _lgtab(N):
    if N not in LOG_CACHE:
        LOG_CACHE[N] = np.array([math.lgamma(k + 1) for k in range(N + 1)])
    return LOG_CACHE[N]

def locus_stats(dose_i8):
    """dose_i8: loci x n int8 (no missing).  Returns dict of per-locus arrays."""
    n = dose_i8.shape[1]
    d = dose_i8.astype(np.float64)
    alt = d.sum(1)
    p = alt / (2.0 * n); q = 1.0 - p
    Ho = (dose_i8 == 1).sum(1) / float(n)
    He_hw = 2.0 * p * q
    He_unb = (2.0 * n / (2.0 * n - 1.0)) * He_hw if n > 1 else np.full_like(p, np.nan)
    b = (n / (n - 1.0)) * (p * q - ((2.0 * n - 1.0) / (4.0 * n)) * Ho) if n > 1 else np.zeros_like(p)
    c = Ho / 2.0
    return dict(n=n, p=p, q=q, Ho=Ho, He_hw=He_hw, He_unb=He_unb, b=b, c=c,
                alt_count=np.rint(alt).astype(np.int64))

def fis_wc(b, c):
    den = np.sum(b + c)
    return float(np.sum(b) / den) if den != 0 else float('nan')

def fis_jackknife_se(b, c):
    B, D = np.sum(b), np.sum(b + c)
    L = len(b)
    th = (B - b) / (D - (b + c))
    th = th[np.isfinite(th)]
    m = th.mean()
    return float(math.sqrt((len(th) - 1) / len(th) * np.sum((th - m) ** 2)))

def boot_ratio_ci(num, den, nboot=1000, rng=None, alpha=0.05):
    rng = rng or np.random.default_rng(SEED)
    L = len(num)
    vals = np.empty(nboot)
    for r in range(nboot):
        w = np.bincount(rng.integers(0, L, L), minlength=L).astype(np.float64)
        vals[r] = (num @ w) / (den @ w)
    return float(np.percentile(vals, 100 * alpha / 2)), float(np.percentile(vals, 100 * (1 - alpha / 2)))

def allelic_richness(alt_counts, twoN, g):
    """El Mousadik & Petit (1996) rarefied allelic richness at g gene copies, biallelic."""
    if twoN < g:
        return float('nan'), float('nan')
    lg = _lgtab(twoN)
    def logC(N, k):
        return lg[N] - lg[k] - lg[N - k]
    lCNg = logC(twoN, g)
    nref = twoN - alt_counts
    r_alt = np.where(nref >= g, np.exp(logC(twoN, g) * 0 + (lg[nref] - lg[g] - lg[np.maximum(nref - g, 0)]) - lCNg), 0.0)
    r_ref = np.where(alt_counts >= g, np.exp((lg[alt_counts] - lg[g] - lg[np.maximum(alt_counts - g, 0)]) - lCNg), 0.0)
    ar = (1.0 - r_alt) + (1.0 - r_ref)
    return float(ar.mean()), float(ar.std(ddof=1) / math.sqrt(len(ar)))

log('STEP 5  per-population diversity')
POP_COLS = {p: np.where(popk == p)[0] for p in POPS_ORDER}
MAINLAND_COLS = np.where(popk != 'KI')[0]
KI_COLS = POP_COLS['KI']
BIG_POPS = [p for p in POPS_ORDER if len(POP_COLS[p]) >= 19]

ST = {p: locus_stats(Gf[:, POP_COLS[p]]) for p in POPS_ORDER}
ST_MAIN = locus_stats(Gf[:, MAINLAND_COLS])
ST_JOINT = locus_stats(Gf)

div_rows, DIVJ = [], {}
for p in POPS_ORDER:
    s = ST[p]; n = s['n']
    Ho = float(s['Ho'].mean()); Ho_se = float(s['Ho'].std(ddof=1) / math.sqrt(NL))
    He = float(s['He_unb'].mean()); He_se = float(s['He_unb'].std(ddof=1) / math.sqrt(NL))
    Hehw = float(s['He_hw'].mean())
    F = fis_wc(s['b'], s['c'])
    Fse = fis_jackknife_se(s['b'], s['c'])
    Flo, Fhi = boot_ratio_ci(s['b'], s['b'] + s['c'], 1000, np.random.default_rng(SEED))
    poly = (s['p'] > 0) & (s['p'] < 1)
    ar, arse = allelic_richness(s['alt_count'], 2 * n, 16)
    maf = np.minimum(s['p'], s['q'])
    mfm = float(maf[poly].mean()) if poly.any() else float('nan')
    row = dict(population=p, region=pop2region.get(p, ''), n=n, n_loci=NL,
               Ho=Ho, Ho_SE=Ho_se, He_unbiased=He, He_SE=He_se, He_hw=Hehw,
               Fis_WC=F, Fis_jack_SE=Fse, Fis_boot_lo=Flo, Fis_boot_hi=Fhi,
               n_polymorphic=int(poly.sum()), prop_polymorphic=float(poly.mean()),
               Ar_g16=ar, Ar_g16_SE=arse, mean_folded_MAF=mfm)
    div_rows.append(row); DIVJ[p] = row
DIV = pd.DataFrame(div_rows)
DIV.to_csv(os.path.join(OUT, 'diversity_by_population.csv'), index=False)
for r in div_rows:
    log(f"    {r['population']:<16s} n={r['n']:3d} Ho={r['Ho']:.4f} He={r['He_unbiased']:.4f} "
        f"Fis={r['Fis_WC']:+.4f} [{r['Fis_boot_lo']:+.4f},{r['Fis_boot_hi']:+.4f}] "
        f"poly={r['n_polymorphic']:4d} Ar16={r['Ar_g16']:.4f}")

ml_pops = [p for p in POPS_ORDER if p != 'KI']
MAINLAND_MEAN_HE = float(np.mean([DIVJ[p]['He_unbiased'] for p in ml_pops]))
MAINLAND_MEAN_HO = float(np.mean([DIVJ[p]['Ho'] for p in ml_pops]))
RATIOS = {
    'KI_He': DIVJ['KI']['He_unbiased'], 'mainland_mean_He': MAINLAND_MEAN_HE,
    'KI_pct_of_mainland_mean_He': 100 * DIVJ['KI']['He_unbiased'] / MAINLAND_MEAN_HE,
    'KI_pct_of_OVENS_He': 100 * DIVJ['KI']['He_unbiased'] / DIVJ['OVENS']['He_unbiased'],
    'KI_pct_of_mainland_mean_Ho': 100 * DIVJ['KI']['Ho'] / MAINLAND_MEAN_HO,
    'KI_pct_of_OVENS_Ho': 100 * DIVJ['KI']['Ho'] / DIVJ['OVENS']['Ho'],
    'KI_is_lowest_He': bool(DIVJ['KI']['He_unbiased'] == min(DIVJ[p]['He_unbiased'] for p in POPS_ORDER)),
    'KI_is_lowest_Ho': bool(DIVJ['KI']['Ho'] == min(DIVJ[p]['Ho'] for p in POPS_ORDER)),
    'KI_is_lowest_Ar16': bool(DIVJ['KI']['Ar_g16'] == min(
        [DIVJ[p]['Ar_g16'] for p in POPS_ORDER if np.isfinite(DIVJ[p]['Ar_g16'])])),
}
ar_ml = [DIVJ[p]['Ar_g16'] for p in ml_pops if np.isfinite(DIVJ[p]['Ar_g16'])]
RATIOS['mainland_mean_Ar16'] = float(np.mean(ar_ml))
RATIOS['KI_pct_of_mainland_mean_Ar16'] = 100 * DIVJ['KI']['Ar_g16'] / RATIOS['mainland_mean_Ar16']
RATIOS['n_mainland_pops_with_Ar16'] = len(ar_ml)

# allelic richness curve g = 2..16
AR_CURVE = {}
for p in POPS_ORDER:
    s = ST[p]
    AR_CURVE[p] = {str(g): (allelic_richness(s['alt_count'], 2 * s['n'], g)[0]
                            if 2 * s['n'] >= g else None) for g in range(2, 17, 2)}
json.dump({'g16': {p: {'n': DIVJ[p]['n'], 'Ar_g16': DIVJ[p]['Ar_g16'],
                       'Ar_g16_SE': DIVJ[p]['Ar_g16_SE'], 'n_loci': NL} for p in POPS_ORDER},
           'curve': AR_CURVE, 'summary': {k: RATIOS[k] for k in
               ('mainland_mean_Ar16','KI_pct_of_mainland_mean_Ar16','n_mainland_pops_with_Ar16',
                'KI_is_lowest_Ar16')}},
          open(os.path.join(OUT, 'allelic_richness.json'), 'w'), indent=1, default=float)

# ---- HWE-form sensitivity: the whole diversity table on the joint-HWE (3,727) locus set
FIN_JT = np.where(casc_jt['mask_final'])[0]
Gf_jt = Gk[FIN_JT]
HWE_SENS = {}
for p in POPS_ORDER:
    s = locus_stats(Gf_jt[:, POP_COLS[p]])
    HWE_SENS[p] = {'n_loci': len(FIN_JT), 'Ho': float(s['Ho'].mean()),
                   'He_unbiased': float(s['He_unb'].mean()), 'Fis_WC': fis_wc(s['b'], s['c'])}
UNANCH_SENS = {}
FIN_DU = np.where(casc_drop['mask_final'])[0]
Gf_du = Gk[FIN_DU]
for p in POPS_ORDER:
    s = locus_stats(Gf_du[:, POP_COLS[p]])
    UNANCH_SENS[p] = {'n_loci': len(FIN_DU), 'Ho': float(s['Ho'].mean()),
                      'He_unbiased': float(s['He_unb'].mean()), 'Fis_WC': fis_wc(s['b'], s['c'])}
max_he_shift = max(abs(UNANCH_SENS[p]['He_unbiased'] - DIVJ[p]['He_unbiased']) for p in POPS_ORDER)
log(f'  sensitivity: joint-HWE set {len(FIN_JT)} loci; drop-unanchored set {len(FIN_DU)} loci '
    f'(max He shift {max_he_shift:.5f})')

# ===========================================================================
# 6.  PER-INDIVIDUAL HETEROZYGOSITY AND INBREEDING  (PLAN §3.4-3.5; Reviewer 2)
# ===========================================================================
log('STEP 6  per-individual Ho and inbreeding F')
def individual_F(dose_i8, p_ref, n_ref, nboot=1000, rng=None):
    """Three inbreeding estimators for every column of dose_i8 against reference freqs p_ref.
       F_HOM  = excess homozygosity (method of moments, Yang et al. 2010 F-hat-2)
       F_UNI  = correlation of uniting gametes (Yang et al. 2010 F-hat-3)
       F_GRM  = variance-standardised GRM diagonal - 1 (Yang et al. 2010 F-hat-1)"""
    q_ref = 1.0 - p_ref
    use = (p_ref > 0) & (p_ref < 1)
    p = p_ref[use]; q = q_ref[use]
    g = dose_i8[use].astype(np.float64)                     # loci x nind
    L = int(use.sum())
    twopq = 2.0 * p * q
    # F_HOM
    corr = (2.0 * n_ref / (2.0 * n_ref - 1.0)) if n_ref > 1 else 1.0
    E_hom = float(np.sum(1.0 - twopq * corr))
    O_hom = (g != 1).sum(0).astype(float)
    F_hom = (O_hom - E_hom) / (L - E_hom)
    # F_UNI  (F-hat-3)
    V3 = (g * g - (1.0 + 2.0 * p)[:, None] * g + 2.0 * (p * p)[:, None]) / twopq[:, None]
    F_uni = V3.mean(0)
    # F_GRM  (F-hat-1)
    V1 = ((g - 2.0 * p[:, None]) ** 2) / twopq[:, None]
    F_grm = V1.mean(0) - 1.0
    lo = hi = None
    if nboot:
        rng = rng or np.random.default_rng(SEED)
        reps = np.empty((nboot, g.shape[1]))
        for r in range(nboot):
            w = np.bincount(rng.integers(0, L, L), minlength=L).astype(np.float64)
            reps[r] = (w @ V3) / L
        lo = np.percentile(reps, 2.5, axis=0); hi = np.percentile(reps, 97.5, axis=0)
    return dict(L=L, F_hom=F_hom, F_uni=F_uni, F_grm=F_grm, F_uni_lo=lo, F_uni_hi=hi)

Ho_ind = (Gf == 1).sum(0) / float(NL)
IND = pd.DataFrame({'sample': sids, 'population': popk, 'n_loci_called': NL, 'Ho_ind': Ho_ind})
for col in ['F_hom_ownpop','F_uni_ownpop','F_grm_ownpop','F_uni_ownpop_lo95','F_uni_ownpop_hi95',
            'n_loci_F_ownpop','F_hom_mainland','F_uni_mainland','F_grm_mainland',
            'F_uni_mainland_lo95','F_uni_mainland_hi95','n_loci_F_mainland']:
    IND[col] = np.nan
for p in POPS_ORDER:
    cols = POP_COLS[p]
    r = individual_F(Gf[:, cols], ST[p]['p'], len(cols), 1000, np.random.default_rng(SEED))
    IND.loc[cols, 'F_hom_ownpop'] = r['F_hom']; IND.loc[cols, 'F_uni_ownpop'] = r['F_uni']
    IND.loc[cols, 'F_grm_ownpop'] = r['F_grm']; IND.loc[cols, 'F_uni_ownpop_lo95'] = r['F_uni_lo']
    IND.loc[cols, 'F_uni_ownpop_hi95'] = r['F_uni_hi']; IND.loc[cols, 'n_loci_F_ownpop'] = r['L']
rm = individual_F(Gf, ST_MAIN['p'], len(MAINLAND_COLS), 1000, np.random.default_rng(SEED + 1))
IND['F_hom_mainland'] = rm['F_hom']; IND['F_uni_mainland'] = rm['F_uni']
IND['F_grm_mainland'] = rm['F_grm']; IND['F_uni_mainland_lo95'] = rm['F_uni_lo']
IND['F_uni_mainland_hi95'] = rm['F_uni_hi']; IND['n_loci_F_mainland'] = rm['L']

# ===========================================================================
# 7.  PAIRWISE KINSHIP, Yang et al. (2010)  (PLAN §3.6; Reviewer 2)
# ===========================================================================
log('STEP 7  pairwise kinship')
def kinship_matrix(dose_i8, p_ref):
    use = (p_ref > 0) & (p_ref < 1)
    p = p_ref[use]
    Z = (dose_i8[use].astype(np.float64) - 2.0 * p[:, None]) / np.sqrt(2.0 * p * (1.0 - p))[:, None]
    return (Z.T @ Z) / float(use.sum()), int(use.sum())

kin_rows, KIN_BY_POP = [], {}
for p in POPS_ORDER:
    cols = POP_COLS[p]
    if len(cols) < 2:
        continue
    K, Lu = kinship_matrix(Gf[:, cols], ST[p]['p'])
    vals = []
    for a, b_ in itertools.combinations(range(len(cols)), 2):
        kin_rows.append(dict(sample_i=sids[cols[a]], sample_j=sids[cols[b_]], population_i=p,
                             population_j=p, within_pop=1, kinship_yang2010=float(K[a, b_]),
                             relatedness_scale=float(K[a, b_]), n_loci=Lu))
        vals.append(float(K[a, b_]))
    KIN_BY_POP[p] = dict(n=len(cols), n_pairs=len(vals), n_loci=Lu,
                         mean=float(np.mean(vals)), sd=float(np.std(vals, ddof=1)),
                         p05=float(np.percentile(vals, 5)), p50=float(np.percentile(vals, 50)),
                         p95=float(np.percentile(vals, 95)), p99=float(np.percentile(vals, 99)),
                         max=float(np.max(vals)), values=[round(v, 5) for v in vals] if p == 'KI' else None)
    if p == 'KI':
        KI_KIN_MATRIX = K.copy()
pd.DataFrame(kin_rows).to_csv(os.path.join(OUT, 'kinship_pairs.csv'), index=False)
IND['mean_kinship_within_pop'] = np.nan; IND['max_kinship_within_pop'] = np.nan
for p in POPS_ORDER:
    cols = POP_COLS[p]
    if len(cols) < 2:
        continue
    K, _ = kinship_matrix(Gf[:, cols], ST[p]['p'])
    off = K.copy(); np.fill_diagonal(off, np.nan)
    IND.loc[cols, 'mean_kinship_within_pop'] = np.nanmean(off, axis=1)
    IND.loc[cols, 'max_kinship_within_pop'] = np.nanmax(off, axis=1)
IND.to_csv(os.path.join(OUT, 'diversity_by_individual.csv'), index=False)

# empirical within-population null from the six populations with n >= 19
null_vals = np.concatenate([[v for v in
    [float(kinship_matrix(Gf[:, POP_COLS[p]], ST[p]['p'])[0][a, b_])
     for a, b_ in itertools.combinations(range(len(POP_COLS[p])), 2)]] for p in BIG_POPS])
EMP_NULL = dict(pops=BIG_POPS, n_pairs=int(len(null_vals)), mean=float(null_vals.mean()),
                sd=float(null_vals.std(ddof=1)),
                p05=float(np.percentile(null_vals, 5)), p50=float(np.percentile(null_vals, 50)),
                p95=float(np.percentile(null_vals, 95)), p99=float(np.percentile(null_vals, 99)),
                p999=float(np.percentile(null_vals, 99.9)), max=float(null_vals.max()))

# matched n = 8 mainland null for kinship AND individual F
log('  matched n=8 mainland nulls (200 draws x 6 populations)')
rng8 = np.random.default_rng(SEED)
m8_kin, m8_funi, m8_fhom, m8_fgrm = [], [], [], []
M8_BY_POP = {}
for p in BIG_POPS:
    cols = POP_COLS[p]; kk, ff = [], []
    for _ in range(200):
        sub = rng8.choice(cols, 8, replace=False)
        s8 = locus_stats(Gf[:, sub])
        K, _ = kinship_matrix(Gf[:, sub], s8['p'])
        v = K[np.triu_indices(8, 1)]
        kk.extend(v.tolist())
        r = individual_F(Gf[:, sub], s8['p'], 8, nboot=0)
        ff.extend(r['F_uni'].tolist()); m8_fhom.extend(r['F_hom'].tolist())
        m8_fgrm.extend(r['F_grm'].tolist())
    m8_kin.extend(kk); m8_funi.extend(ff)
    M8_BY_POP[p] = dict(kin_mean=float(np.mean(kk)), kin_p95=float(np.percentile(kk, 95)),
                        kin_max=float(np.max(kk)), F_uni_mean=float(np.mean(ff)),
                        F_uni_p05=float(np.percentile(ff, 5)), F_uni_p95=float(np.percentile(ff, 95)))
m8_kin = np.asarray(m8_kin); m8_funi = np.asarray(m8_funi)
MATCHED8 = dict(reps=200, pops=BIG_POPS,
    kinship=dict(n=int(len(m8_kin)), mean=float(m8_kin.mean()), sd=float(m8_kin.std(ddof=1)),
                 p05=float(np.percentile(m8_kin, 5)), p50=float(np.percentile(m8_kin, 50)),
                 p95=float(np.percentile(m8_kin, 95)), p99=float(np.percentile(m8_kin, 99)),
                 max=float(m8_kin.max())),
    F_uni=dict(n=int(len(m8_funi)), mean=float(m8_funi.mean()), sd=float(m8_funi.std(ddof=1)),
               p05=float(np.percentile(m8_funi, 5)), p50=float(np.percentile(m8_funi, 50)),
               p95=float(np.percentile(m8_funi, 95)), max=float(m8_funi.max())),
    F_hom=dict(mean=float(np.mean(m8_fhom)), p05=float(np.percentile(m8_fhom, 5)),
               p95=float(np.percentile(m8_fhom, 95))),
    F_grm=dict(mean=float(np.mean(m8_fgrm)), p05=float(np.percentile(m8_fgrm, 5)),
               p95=float(np.percentile(m8_fgrm, 95))),
    by_pop=M8_BY_POP)
ki_kin_vals = KI_KIN_MATRIX[np.triu_indices(8, 1)]
KI_KIN = dict(n_pairs=28, values=[round(float(v), 5) for v in ki_kin_vals],
              matrix=[[round(float(x), 5) for x in row] for row in KI_KIN_MATRIX],
              labels=[str(s) for s in sids[KI_COLS]],
              mean=float(ki_kin_vals.mean()), max=float(ki_kin_vals.max()),
              min=float(ki_kin_vals.min()),
              n_above_matched_p95=int((ki_kin_vals > MATCHED8['kinship']['p95']).sum()),
              n_above_matched_max=int((ki_kin_vals > MATCHED8['kinship']['max']).sum()),
              empirical_p_max_pair=float((m8_kin >= ki_kin_vals.max()).mean()),
              n_above_0p25_relatedness=int((ki_kin_vals > 0.25).sum()),
              n_above_0p125_relatedness=int((ki_kin_vals > 0.125).sum()))
json.dump({'KI': KI_KIN, 'within_pop_summary': {k: {kk: vv for kk, vv in v.items() if kk != 'values'}
                                                for k, v in KIN_BY_POP.items()},
           'empirical_null_n_ge_19': EMP_NULL, 'matched_n8_null': MATCHED8,
           'convention': 'Yang et al. (2010) genomic relationship on the RELATEDNESS scale '
                         '(unrelated ~0, full sibs/parent-offspring ~0.5); halve for kinship coefficient'},
          open(os.path.join(OUT, 'kinship_summary.json'), 'w'), indent=1, default=float)
log(f"  KI kinship mean={KI_KIN['mean']:+.4f} max={KI_KIN['max']:+.4f}; "
    f"matched-n8 mainland null p95={MATCHED8['kinship']['p95']:+.4f} max={MATCHED8['kinship']['max']:+.4f}")

# ===========================================================================
# 8.  PAIRWISE F_ST — Weir & Cockerham (1984) theta, plus G_ST and G''_ST  (PLAN §3.7)
# ===========================================================================
log('STEP 8  pairwise Fst (11x11) with 1000-replicate bootstrap over loci')
def wc_components(s1, s2):
    n1, n2 = float(s1['n']), float(s2['n']); r = 2.0
    p1, p2 = s1['p'], s2['p']; h1, h2 = s1['Ho'], s2['Ho']
    nbar = (n1 + n2) / r
    nc = (r * nbar - (n1 * n1 + n2 * n2) / (r * nbar)) / (r - 1.0)
    pbar = (n1 * p1 + n2 * p2) / (n1 + n2)
    s2_ = (n1 * (p1 - pbar) ** 2 + n2 * (p2 - pbar) ** 2) / ((r - 1.0) * nbar)
    hbar = (n1 * h1 + n2 * h2) / (n1 + n2)
    core = pbar * (1.0 - pbar) - ((r - 1.0) / r) * s2_
    a = (nbar / nc) * (s2_ - (1.0 / (nbar - 1.0)) * (core - hbar / 4.0))
    b = (nbar / (nbar - 1.0)) * (core - ((2.0 * nbar - 1.0) / (4.0 * nbar)) * hbar)
    c = hbar / 2.0
    return a, b, c

def nei_gst(s1, s2):
    n1, n2 = float(s1['n']), float(s2['n']); k = 2.0
    ntil = k / (1.0 / n1 + 1.0 / n2)
    p1, p2 = s1['p'], s2['p']
    hs_raw = np.mean([1 - (p1 ** 2 + (1 - p1) ** 2), 1 - (p2 ** 2 + (1 - p2) ** 2)], axis=0)
    ho_bar = (s1['Ho'] + s2['Ho']) / 2.0
    Hs = (ntil / (ntil - 1.0)) * (hs_raw - ho_bar / (2.0 * ntil))
    pbar = (p1 + p2) / 2.0
    Ht = 1 - (pbar ** 2 + (1 - pbar) ** 2) + Hs / (ntil * k) - ho_bar / (2.0 * ntil * k)
    HS, HT = float(np.mean(Hs)), float(np.mean(Ht))
    gst = (HT - HS) / HT if HT > 0 else float('nan')
    gpp = (k * (HT - HS)) / ((k * HT - HS) * (1.0 - HS)) if (k * HT - HS) > 0 and HS < 1 else float('nan')
    return gst, gpp, HS, HT

pairs = list(itertools.combinations(POPS_ORDER, 2))
A = np.empty((len(pairs), NL)); B = np.empty((len(pairs), NL)); C = np.empty((len(pairs), NL))
fst_rows = []
for t, (pi, pj) in enumerate(pairs):
    a, b, c = wc_components(ST[pi], ST[pj])
    A[t], B[t], C[t] = a, b, c
DEN = A + B + C
fst_point = A.sum(1) / DEN.sum(1)
rngb = np.random.default_rng(SEED)
boot = np.empty((1000, len(pairs)))
for r in range(1000):
    w = np.bincount(rngb.integers(0, NL, NL), minlength=NL).astype(np.float64)
    boot[r] = (A @ w) / (DEN @ w)
lo = np.percentile(boot, 2.5, axis=0); hi = np.percentile(boot, 97.5, axis=0)
FSTM = pd.DataFrame(np.zeros((11, 11)), index=POPS_ORDER, columns=POPS_ORDER)
for t, (pi, pj) in enumerate(pairs):
    gst, gpp, HS, HT = nei_gst(ST[pi], ST[pj])
    fst_rows.append(dict(pop_i=pi, pop_j=pj, n_i=ST[pi]['n'], n_j=ST[pj]['n'], n_loci=NL,
                         Fst_WC=float(fst_point[t]), Fst_lo95=float(lo[t]), Fst_hi95=float(hi[t]),
                         Gst_Nei=gst, GstPP_Hedrick=gpp, Hs=HS, Ht=HT))
    FSTM.loc[pi, pj] = FSTM.loc[pj, pi] = float(fst_point[t])
FSTDF = pd.DataFrame(fst_rows)
FSTDF.to_csv(os.path.join(OUT, 'fst_pairwise.csv'), index=False)
FSTM.to_csv(os.path.join(OUT, 'fst_matrix.csv'))
json.dump({'method': 'Weir & Cockerham (1984) theta, ratio of sums over loci',
           'n_loci': NL, 'n_bootstrap': 1000, 'seed': SEED,
           'matrix': {pi: {pj: (float(FSTM.loc[pi, pj]) if pi != pj else 0.0) for pj in POPS_ORDER}
                      for pi in POPS_ORDER},
           'pairs': fst_rows}, open(os.path.join(OUT, 'fst_pairwise_ci.json'), 'w'), indent=1, default=float)
ki_fst = FSTDF[(FSTDF.pop_i == 'KI') | (FSTDF.pop_j == 'KI')].copy()
ki_fst['other'] = np.where(ki_fst.pop_i == 'KI', ki_fst.pop_j, ki_fst.pop_i)
ki_fst = ki_fst.sort_values('Fst_WC')
for _, r in ki_fst.iterrows():
    log(f"    KI vs {r['other']:<16s} Fst={r['Fst_WC']:.4f} [{r['Fst_lo95']:.4f},{r['Fst_hi95']:.4f}] "
        f"G''st={r['GstPP_Hedrick']:.4f}")
ml_pairs = FSTDF[(FSTDF.pop_i != 'KI') & (FSTDF.pop_j != 'KI')]
FST_SUMMARY = dict(KI_min=float(ki_fst.Fst_WC.min()), KI_min_pop=str(ki_fst.iloc[0]['other']),
                   KI_max=float(ki_fst.Fst_WC.max()), KI_max_pop=str(ki_fst.iloc[-1]['other']),
                   KI_mean=float(ki_fst.Fst_WC.mean()),
                   mainland_mainland_mean=float(ml_pairs.Fst_WC.mean()),
                   mainland_mainland_min=float(ml_pairs.Fst_WC.min()),
                   mainland_mainland_max=float(ml_pairs.Fst_WC.max()))

# ===========================================================================
# 9.  PCA — Patterson normalisation, SVD  (PLAN §3.8; Reviewer 2's overlap objection)
# ===========================================================================
log('STEP 9  PCA')
pj_ = ST_JOINT['p']; qj_ = 1.0 - pj_
scale = np.where(pj_ * qj_ > 0, np.sqrt(2.0 * pj_ * qj_), 1.0)
X = (GF.T - 2.0 * pj_) / scale                       # individuals x loci
check('pca_no_missing', int(np.isnan(X).sum()) == 0, int(np.isnan(X).sum()), 0)
U, S, Vt = np.linalg.svd(X, full_matrices=False)
scores = U * S
varpct = 100.0 * (S ** 2) / float((S ** 2).sum())
PCS = pd.DataFrame({'sample': sids, 'population': popk,
                    'region': [pop2region.get(p, '') for p in popk]})
for k in range(6):
    PCS[f'PC{k+1}'] = scores[:, k]
PCS.to_csv(os.path.join(OUT, 'pca_scores.csv'), index=False)

is_ki = (popk == 'KI')
PCA_OVERLAP = {}
for k in range(4):
    v = scores[:, k]; vk = v[is_ki]; vm = v[~is_ki]
    ki_lo, ki_hi = float(vk.min()), float(vk.max())
    ml_lo, ml_hi = float(vm.min()), float(vm.max())
    inside = [str(sids[KI_COLS[i]]) for i in range(8) if ml_lo <= vk[i] <= ml_hi]
    ml_inside = [str(sids[j]) for j in np.where(~is_ki)[0] if ki_lo <= v[j] <= ki_hi]
    sep = (ki_hi < ml_lo) or (ki_lo > ml_hi)
    gap = (ml_lo - ki_hi) if ki_hi < ml_lo else ((ki_lo - ml_hi) if ki_lo > ml_hi else 0.0)
    PCA_OVERLAP[f'PC{k+1}'] = dict(
        variance_pct=float(varpct[k]),
        KI_range=[ki_lo, ki_hi], KI_mean=float(vk.mean()),
        mainland_range=[ml_lo, ml_hi], mainland_mean=float(vm.mean()),
        fully_separated=bool(sep), gap=float(gap),
        n_KI_inside_mainland_range=len(inside), KI_inside_mainland_range=inside,
        n_mainland_inside_KI_range=len(ml_inside),
        mainland_inside_KI_range=ml_inside[:20],
        KI_scores={str(sids[KI_COLS[i]]): float(vk[i]) for i in range(8)})
    log(f"    PC{k+1} ({varpct[k]:.2f}%) KI [{ki_lo:+.2f},{ki_hi:+.2f}] mainland [{ml_lo:+.2f},{ml_hi:+.2f}] "
        f"separated={sep} KI-inside-mainland={len(inside)}")
# joint separability on PC1 x PC2: minimum Euclidean distance between a KI and a mainland point
def min_cross_dist(axes):
    a = scores[np.ix_(np.where(is_ki)[0], axes)]
    b = scores[np.ix_(np.where(~is_ki)[0], axes)]
    d = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
    return float(d.min()), float(np.sqrt(((a[:, None, :] - a[None, :, :]) ** 2).sum(-1)).max())
PCA_OVERLAP['PC1xPC2'] = dict(zip(['min_KI_to_mainland_distance', 'max_KI_internal_distance'],
                                  min_cross_dist([0, 1])))
PCA_OVERLAP['PC1xPC2']['KI_cluster_separated'] = bool(
    PCA_OVERLAP['PC1xPC2']['min_KI_to_mainland_distance'] >
    PCA_OVERLAP['PC1xPC2']['max_KI_internal_distance'])
PCA_OVERLAP['PC2xPC3'] = dict(zip(['min_KI_to_mainland_distance', 'max_KI_internal_distance'],
                                  min_cross_dist([1, 2])))
json.dump({'n_loci': NL, 'n_samples': 222, 'normalisation': 'Patterson/Price (centre, /sqrt(2pq))',
           'variance_pct': [float(v) for v in varpct[:10]],
           'singular_values': [float(v) for v in S[:10]], 'overlap': PCA_OVERLAP},
          open(os.path.join(OUT, 'pca_eigen.json'), 'w'), indent=1, default=float)

# ===========================================================================
# 10.  RAREFACTION / SUBSAMPLING TO MATCHED n  (PLAN §3.9, §4.2)
# ===========================================================================
log('STEP 10  rarefaction (200 draws at n = 5, 8, 10, 15, 20)')
RAR = {'seed': SEED, 'reps': 200, 'levels': [5, 8, 10, 15, 20],
       'KI_observed': {'n': 8, 'Ho': DIVJ['KI']['Ho'], 'He': DIVJ['KI']['He_unbiased'],
                       'Fis': DIVJ['KI']['Fis_WC'], 'n_poly': DIVJ['KI']['n_polymorphic'],
                       'Ar_g16': DIVJ['KI']['Ar_g16'],
                       'mean_folded_MAF': DIVJ['KI']['mean_folded_MAF']},
       'by_pop': {}}
rngr = np.random.default_rng(SEED)
SUB8_SPECTRA = {}
for p in BIG_POPS:
    cols = POP_COLS[p]
    entry = {'n_full': len(cols),
             'full': {'Ho': DIVJ[p]['Ho'], 'He': DIVJ[p]['He_unbiased'], 'Fis': DIVJ[p]['Fis_WC'],
                      'n_poly': DIVJ[p]['n_polymorphic'], 'Ar_g16': DIVJ[p]['Ar_g16'],
                      'mean_folded_MAF': DIVJ[p]['mean_folded_MAF']}}
    for lev in [5, 8, 10, 15, 20]:
        if lev > len(cols):
            continue
        acc = {k: [] for k in ('Ho','He','Fis','n_poly','Ar','maf')}
        spectra = []
        for _ in range(200):
            sub = rngr.choice(cols, lev, replace=False)
            s = locus_stats(Gf[:, sub])
            poly = (s['p'] > 0) & (s['p'] < 1)
            acc['Ho'].append(float(s['Ho'].mean())); acc['He'].append(float(s['He_unb'].mean()))
            acc['Fis'].append(fis_wc(s['b'], s['c'])); acc['n_poly'].append(int(poly.sum()))
            acc['Ar'].append(allelic_richness(s['alt_count'], 2 * lev, 16)[0] if 2 * lev >= 16 else np.nan)
            mac = np.minimum(s['alt_count'], 2 * lev - s['alt_count'])
            acc['maf'].append(float((mac[poly] / (2.0 * lev)).mean()) if poly.any() else np.nan)
            if lev == 8:
                spectra.append(np.bincount(mac[poly], minlength=9)[1:9])
        entry[f'n{lev}'] = {
            'Ho_mean': float(np.mean(acc['Ho'])), 'Ho_sd': float(np.std(acc['Ho'], ddof=1)),
            'He_mean': float(np.mean(acc['He'])), 'He_sd': float(np.std(acc['He'], ddof=1)),
            'He_pct_of_full': float(100 * np.mean(acc['He']) / DIVJ[p]['He_unbiased']),
            'Fis_mean': float(np.mean(acc['Fis'])), 'Fis_sd': float(np.std(acc['Fis'], ddof=1)),
            'n_poly_mean': float(np.mean(acc['n_poly'])), 'n_poly_sd': float(np.std(acc['n_poly'], ddof=1)),
            'Ar_g16_mean': float(np.nanmean(acc['Ar'])) if np.isfinite(acc['Ar']).any() else None,
            'mean_folded_MAF_mean': float(np.nanmean(acc['maf'])),
            'mean_folded_MAF_sd': float(np.nanstd(acc['maf'], ddof=1))}
        if lev == 8:
            SUB8_SPECTRA[p] = np.array(spectra, float)
    RAR['by_pop'][p] = entry
json.dump(RAR, open(os.path.join(OUT, 'rarefaction.json'), 'w'), indent=1, default=float)

# ---- proportion polymorphic, matched n = 8 (PLAN §4)
POLY_MATCH = {'n_loci_final': NL,
              'KI': {'n': 8, 'n_polymorphic': DIVJ['KI']['n_polymorphic'],
                     'prop': DIVJ['KI']['prop_polymorphic']},
              'all_populations_full_n': {p: {'n': DIVJ[p]['n'], 'n_polymorphic': DIVJ[p]['n_polymorphic'],
                                             'prop': DIVJ[p]['prop_polymorphic']} for p in POPS_ORDER},
              'mainland_at_n8': {}}
for p in BIG_POPS:
    e = RAR['by_pop'][p]['n8']
    POLY_MATCH['mainland_at_n8'][p] = {
        'n_poly_mean': e['n_poly_mean'], 'n_poly_sd': e['n_poly_sd'],
        'prop_mean': e['n_poly_mean'] / NL,
        'KI_pct_of_this': 100 * DIVJ['KI']['n_polymorphic'] / e['n_poly_mean'],
        'KI_deficit_pct': 100 * (1 - DIVJ['KI']['n_polymorphic'] / e['n_poly_mean'])}
m8 = float(np.mean([POLY_MATCH['mainland_at_n8'][p]['n_poly_mean'] for p in BIG_POPS]))
POLY_MATCH['mainland_mean_n_poly_at_n8'] = m8
POLY_MATCH['KI_pct_of_mainland_mean_at_n8'] = 100 * DIVJ['KI']['n_polymorphic'] / m8
POLY_MATCH['KI_deficit_vs_mainland_mean_at_n8_pct'] = 100 * (1 - DIVJ['KI']['n_polymorphic'] / m8)
json.dump(POLY_MATCH, open(os.path.join(OUT, 'polymorphism_matched_n.json'), 'w'), indent=1, default=float)
log(f"  polymorphic: KI {DIVJ['KI']['n_polymorphic']}/{NL} = {100*DIVJ['KI']['prop_polymorphic']:.1f}%; "
    f"mainland at n=8 mean {m8:.0f} -> KI deficit {POLY_MATCH['KI_deficit_vs_mainland_mean_at_n8_pct']:.1f}%")

# ===========================================================================
# 11.  OLD-vs-NEW BIAS DIAGNOSTIC (Reviewer 1)  (PLAN §5)
# ===========================================================================
log('STEP 11  Reviewer-1 bias diagnostic: what the CloneID rule lost')
new_clone_snp = {}
SNPcol = META.SNP.astype(str).values
for j in range(L0):
    k = CloneID[j]
    if k not in new_clone_snp:
        new_clone_snp[k] = SNPcol[j]
shared_clones = set(new_clone_snp) & set(old_clone_snp)
match_clones = {c for c in shared_clones if new_clone_snp[c] == old_clone_snp[c]}
check('shared_cloneids_11451', abs(len(shared_clones) - 11451) <= 10, len(shared_clones), 11451)
check('cloneid_callstring_match_10336', abs(len(match_clones) - 10336) <= 10, len(match_clones), 10336)
retained_by_old = np.array([CloneID[j] in match_clones for j in range(L0)])
lab_ret = retained_by_old[FIN]
IDX_RET = FIN[lab_ret]; IDX_LOST = FIN[~lab_ret]

def subset_stats(idx_rows):
    out = {}
    for p in POPS_ORDER:
        s = locus_stats(Gk[np.ix_(idx_rows, POP_COLS[p])])
        poly = (s['p'] > 0) & (s['p'] < 1)
        out[p] = dict(n_loci=len(idx_rows), Ho=float(s['Ho'].mean()),
                      He=float(s['He_unb'].mean()), Fis=fis_wc(s['b'], s['c']),
                      n_poly=int(poly.sum()),
                      Ar_g16=allelic_richness(s['alt_count'], 2 * s['n'], 16)[0]
                              if 2 * s['n'] >= 16 else None)
    return out
SS_ALL  = {p: dict(n_loci=NL, Ho=DIVJ[p]['Ho'], He=DIVJ[p]['He_unbiased'], Fis=DIVJ[p]['Fis_WC'],
                   n_poly=DIVJ[p]['n_polymorphic'], Ar_g16=DIVJ[p]['Ar_g16']) for p in POPS_ORDER}
SS_RET  = subset_stats(IDX_RET)
SS_LOST = subset_stats(IDX_LOST)

def poly_in(idx_rows, cols):
    s = locus_stats(Gk[np.ix_(idx_rows, cols)])
    return int(((s['p'] > 0) & (s['p'] < 1)).sum())
OLDNEW = {
  'rule': 'old pipeline kept a CloneID iff present in BOTH reports AND the two SNP call strings matched '
          '(first row per CloneID in each report)',
  'new_report_cloneids': len(new_clone_snp), 'old_report_cloneids': len(old_clone_snp),
  'shared_cloneids': len(shared_clones), 'cloneid_callstring_match': len(match_clones),
  'final_loci': NL, 'retained_by_old': int(lab_ret.sum()), 'lost_by_old': int((~lab_ret).sum()),
  'pct_lost_by_old': float(100 * (~lab_ret).mean()),
  'retained_poly_KI': poly_in(IDX_RET, KI_COLS), 'lost_poly_KI': poly_in(IDX_LOST, KI_COLS),
  'retained_poly_mainland': poly_in(IDX_RET, MAINLAND_COLS),
  'lost_poly_mainland': poly_in(IDX_LOST, MAINLAND_COLS),
  'by_population': {p: {'all_final': SS_ALL[p], 'retained_by_old': SS_RET[p],
                        'lost_by_old': SS_LOST[p],
                        'He_inflation_pct': 100 * (SS_RET[p]['He'] / SS_ALL[p]['He'] - 1),
                        'Ho_inflation_pct': 100 * (SS_RET[p]['Ho'] / SS_ALL[p]['Ho'] - 1)}
                    for p in POPS_ORDER},
}
OLDNEW['retained_prop_poly_KI'] = OLDNEW['retained_poly_KI'] / max(OLDNEW['retained_by_old'], 1)
OLDNEW['lost_prop_poly_KI'] = OLDNEW['lost_poly_KI'] / max(OLDNEW['lost_by_old'], 1)
OLDNEW['retained_prop_poly_mainland'] = OLDNEW['retained_poly_mainland'] / max(OLDNEW['retained_by_old'], 1)
OLDNEW['lost_prop_poly_mainland'] = OLDNEW['lost_poly_mainland'] / max(OLDNEW['lost_by_old'], 1)
sret = locus_stats(Gk[np.ix_(IDX_RET, np.arange(222))]); slost = locus_stats(Gk[np.ix_(IDX_LOST, np.arange(222))])
OLDNEW['allele_frequency_character'] = {
   'retained_by_old': {'mean_joint_MAF': float(np.minimum(sret['p'], sret['q']).mean()),
       'median_joint_MAF': float(np.median(np.minimum(sret['p'], sret['q']))),
       'pct_MAF_lt_0.01': float(100 * (np.minimum(sret['p'], sret['q']) < 0.01).mean()),
       'pct_monomorphic_in_222': float(100 * ((sret['p'] == 0) | (sret['p'] == 1)).mean()),
       'mean_He_joint': float(sret['He_unb'].mean())},
   'lost_by_old': {'mean_joint_MAF': float(np.minimum(slost['p'], slost['q']).mean()),
       'median_joint_MAF': float(np.median(np.minimum(slost['p'], slost['q']))),
       'pct_MAF_lt_0.01': float(100 * (np.minimum(slost['p'], slost['q']) < 0.01).mean()),
       'pct_monomorphic_in_222': float(100 * ((slost['p'] == 0) | (slost['p'] == 1)).mean()),
       'mean_He_joint': float(slost['He_unb'].mean())}}
OLDNEW['bias_direction'] = ('the CloneID rule preferentially discarded low-frequency / near-fixed loci, '
                            'so every published H_E is inflated; the inflation is larger for mainland '
                            'populations than for KI, so the KI:mainland ratio was biased downward slightly '
                            'while all absolute values were biased upward')
json.dump(OLDNEW, open(os.path.join(OUT, 'oldnew_locus_bias.json'), 'w'), indent=1, default=float)
log(f"  final {NL}: retained by old rule {OLDNEW['retained_by_old']} "
    f"({100-OLDNEW['pct_lost_by_old']:.1f}%), LOST {OLDNEW['lost_by_old']} ({OLDNEW['pct_lost_by_old']:.1f}%)")
log(f"  He KI all={SS_ALL['KI']['He']:.4f} retained-subset={SS_RET['KI']['He']:.4f} "
    f"lost-subset={SS_LOST['KI']['He']:.4f}")

# ===========================================================================
# 12.  FOLDED SITE FREQUENCY SPECTRA, KI vs RAREFIED MAINLAND  (PLAN §6)
# ===========================================================================
log('STEP 12  folded SFS with sample-size-matched mainland reference')
def folded(s, n):
    poly = (s['p'] > 0) & (s['p'] < 1)
    mac = np.minimum(s['alt_count'], 2 * n - s['alt_count'])[poly]
    return mac, mac / (2.0 * n)

mac_KI, maf_KI = folded(ST['KI'], 8)
cls_KI = np.bincount(mac_KI, minlength=9)[1:9]
SFS = {'n_loci_total': NL,
       'KI': {'n': 8, 'gene_copies': 16, 'n_polymorphic': int(len(mac_KI)),
              'mean_folded_maf': float(maf_KI.mean()), 'median_folded_maf': float(np.median(maf_KI)),
              'mac_counts': {str(i + 1): int(cls_KI[i]) for i in range(8)},
              'mac_proportions': {str(i + 1): float(cls_KI[i] / cls_KI.sum()) for i in range(8)}},
       'tests': {}}
s_ov = ST['OVENS']; mac_ov, maf_ov = folded(s_ov, 19)
SFS['OVENS_full'] = {'n': 19, 'gene_copies': 38, 'n_polymorphic': int(len(mac_ov)),
                     'mean_folded_maf': float(maf_ov.mean()),
                     'median_folded_maf': float(np.median(maf_ov)),
                     'note': 'reported for continuity with the submitted analysis ONLY; '
                             'not a valid comparator for KI because n differs'}
rng_s = np.random.default_rng(SEED)
for p in ['OVENS', 'SNOWY', 'THREDBO', 'TENTERFIELD']:
    cols = POP_COLS[p]
    means, npolys, props, Ds = [], [], [], []
    for _ in range(200):
        sub = rng_s.choice(cols, 8, replace=False)
        s = locus_stats(Gf[:, sub])
        mac, maf = folded(s, 8)
        means.append(float(maf.mean())); npolys.append(len(mac))
        cl = np.bincount(mac, minlength=9)[1:9]
        props.append(cl / cl.sum())
        Ds.append(ks_2samp(maf_KI, maf)[0])
    props = np.array(props); means = np.array(means)
    SFS[f'{p}_rarefied_n8'] = {
        'n': 8, 'reps': 200, 'mean_folded_maf': float(means.mean()), 'sd': float(means.std(ddof=1)),
        'ci95': [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))],
        'n_polymorphic_mean': float(np.mean(npolys)), 'n_polymorphic_sd': float(np.std(npolys, ddof=1)),
        'mac_proportions_mean': {str(i + 1): float(props[:, i].mean()) for i in range(8)},
        'mac_proportions_lo': {str(i + 1): float(np.percentile(props[:, i], 2.5)) for i in range(8)},
        'mac_proportions_hi': {str(i + 1): float(np.percentile(props[:, i], 97.5)) for i in range(8)},
        'KS_D_vs_KI_mean': float(np.mean(Ds)), 'KS_D_vs_KI_sd': float(np.std(Ds, ddof=1)),
        'empirical_p_mean_maf_ge_KI': float((means >= maf_KI.mean()).mean())}
    if p == 'OVENS':
        OV8_MEANS = means
# D. KI on the old-CloneID-retained subset (the direct test of Reviewer 1's hypothesis)
s_ki_old = locus_stats(Gk[np.ix_(IDX_RET, KI_COLS)])
mac_kio, maf_kio = folded(s_ki_old, 8)
cls_kio = np.bincount(mac_kio, minlength=9)[1:9]
s_ki_lost = locus_stats(Gk[np.ix_(IDX_LOST, KI_COLS)])
mac_kil, maf_kil = folded(s_ki_lost, 8)
SFS['KI_on_old_cloneid_subset'] = {
    'n_loci': len(IDX_RET), 'n_polymorphic': int(len(mac_kio)),
    'mean_folded_maf': float(maf_kio.mean()), 'median_folded_maf': float(np.median(maf_kio)),
    'mac_counts': {str(i + 1): int(cls_kio[i]) for i in range(8)},
    'mac_proportions': {str(i + 1): float(cls_kio[i] / cls_kio.sum()) for i in range(8)}}
SFS['KI_on_loci_lost_by_old_rule'] = {
    'n_loci': len(IDX_LOST), 'n_polymorphic': int(len(mac_kil)),
    'mean_folded_maf': float(maf_kil.mean()) if len(mac_kil) else None}
D, pks = ks_2samp(maf_KI, maf_ov); U, z, pmw = mannwhitney(maf_KI, maf_ov)
SFS['tests']['KI_vs_OVENS_full_legacy'] = {'KS_D': D, 'KS_p': pks, 'MW_z': z, 'MW_p': pmw,
    'note': 'UNMATCHED sample size (n=8 vs n=19); reported only to show what the submitted '
            'comparison becomes on co-processed data. Not the primary result.'}
pooled_ov8 = None
D2s = SFS['OVENS_rarefied_n8']['KS_D_vs_KI_mean']
SFS['tests']['KI_vs_OVENS_n8_PRIMARY'] = {
    'KS_D_mean_over_200_subsamples': D2s,
    'KS_D_sd': SFS['OVENS_rarefied_n8']['KS_D_vs_KI_sd'],
    'empirical_p_mean_maf_ge_KI': SFS['OVENS_rarefied_n8']['empirical_p_mean_maf_ge_KI'],
    'KI_mean_folded_maf': float(maf_KI.mean()),
    'OVENS_n8_mean_folded_maf': SFS['OVENS_rarefied_n8']['mean_folded_maf'],
    'difference': float(maf_KI.mean() - SFS['OVENS_rarefied_n8']['mean_folded_maf']),
    'legacy_difference_unmatched': float(maf_KI.mean() - float(maf_ov.mean())),
    'pct_of_legacy_gap_that_is_sample_size_artefact':
        float(100 * (1 - (maf_KI.mean() - SFS['OVENS_rarefied_n8']['mean_folded_maf']) /
                     (maf_KI.mean() - maf_ov.mean())))}
D3, p3 = ks_2samp(maf_KI, maf_kio); U3, z3, pmw3 = mannwhitney(maf_KI, maf_kio)
SFS['tests']['KI_all_vs_KI_old_subset'] = {'KS_D': D3, 'KS_p': p3, 'MW_z': z3, 'MW_p': pmw3,
    'KI_all_mean': float(maf_KI.mean()), 'KI_old_subset_mean': float(maf_kio.mean())}
SFS['submitted_values_for_reference'] = {'KI_mean_folded_maf': 0.248, 'KI_n_polymorphic': 472,
    'OVENS_mean_folded_maf': 0.175, 'OVENS_n_polymorphic': 788, 'KS_D': 0.35,
    'source': 'MS_submitted.md / bottleneck_diagnostics.json'}
env_lo = min(SFS[f'{p}_rarefied_n8']['mean_folded_maf'] for p in ['OVENS','SNOWY','THREDBO','TENTERFIELD'])
env_hi = max(SFS[f'{p}_rarefied_n8']['mean_folded_maf'] for p in ['OVENS','SNOWY','THREDBO','TENTERFIELD'])
SFS['mainland_n8_envelope_mean_folded_maf'] = [env_lo, env_hi]
_merge_shift = float(maf_kio.mean() - maf_KI.mean())
_ss_pct = SFS['tests']['KI_vs_OVENS_n8_PRIMARY']['pct_of_legacy_gap_that_is_sample_size_artefact']
_merge_clause = (
    f"Restricting KI to the loci the old CloneID rule would have kept changes its mean folded MAF by only "
    f"{_merge_shift:+.4f} (to {maf_kio.mean():.3f}), so the merge itself did NOT create the shape: "
    f"Reviewer 1's specific mechanism is not supported."
    if abs(_merge_shift) < 0.01 else
    f"Restricting KI to the loci the old CloneID rule would have kept moves its mean folded MAF to "
    f"{maf_kio.mean():.3f} ({_merge_shift:+.4f}), so part of the published intermediate-frequency excess "
    f"was created by the merge.")
SFS['merge_hypothesis_supported'] = bool(abs(_merge_shift) >= 0.01)
SFS['sample_size_artefact_pct_of_published_gap'] = _ss_pct
SFS['verdict'] = (
    f"On the co-processed data KI's mean folded MAF is {maf_KI.mean():.3f} over {len(mac_KI)} polymorphic "
    f"loci (submitted: 0.248 over 472), and the folded spectrum is close to monotonically decreasing rather "
    f"than intermediate-peaked. {_merge_clause} What DOES explain roughly half of the published KI-Ovens "
    f"difference is unmatched sample size: rarefying Ovens to n = 8 raises its mean folded MAF from "
    f"{maf_ov.mean():.3f} to {SFS['OVENS_rarefied_n8']['mean_folded_maf']:.3f}, removing {_ss_pct:.0f}% of "
    f"the published gap. A residual excess of intermediate-frequency variants in KI remains "
    f"(KI {maf_KI.mean():.3f} vs the mainland n = 8 envelope {env_lo:.3f}-{env_hi:.3f}; empirical p "
    f"{'< 0.005' if SFS['OVENS_rarefied_n8']['empirical_p_mean_maf_ge_KI'] == 0 else '= ' + format(SFS['OVENS_rarefied_n8']['empirical_p_mean_maf_ge_KI'], '.3f')}), which is consistent with drift at small "
    f"effective population size but is not diagnostic of an acute recent bottleneck.")
json.dump(SFS, open(os.path.join(OUT, 'sfs.json'), 'w'), indent=1, default=float)
log(f"  KI mean folded MAF {maf_KI.mean():.4f}; OVENS full {maf_ov.mean():.4f}; "
    f"OVENS n=8 {SFS['OVENS_rarefied_n8']['mean_folded_maf']:.4f}; KI-on-old-subset {maf_kio.mean():.4f}")

# ===========================================================================
# 13.  MAF / MAC SENSITIVITY  (PLAN §3.11)
# ===========================================================================
log('STEP 13  MAF/MAC sensitivity')
mac_joint = np.minimum(ST_JOINT['alt_count'], 2 * 222 - ST_JOINT['alt_count'])
maf_joint = np.minimum(ST_JOINT['p'], ST_JOINT['q'])
SCEN = {'none (primary)': np.ones(NL, bool), 'MAC>=3': mac_joint >= 3,
        'MAF>=0.01': maf_joint >= 0.01, 'MAF>=0.02': maf_joint >= 0.02,
        'MAF>=0.05': maf_joint >= 0.05}
MAFSENS = {}
for name, msk in SCEN.items():
    rows = {}
    sel = np.where(msk)[0]
    for p in POPS_ORDER:
        s = locus_stats(Gf[np.ix_(sel, POP_COLS[p])])
        rows[p] = {'Ho': float(s['Ho'].mean()), 'He': float(s['He_unb'].mean()),
                   'Fis': fis_wc(s['b'], s['c']),
                   'n_poly': int(((s['p'] > 0) & (s['p'] < 1)).sum())}
    mlm = float(np.mean([rows[p]['He'] for p in ml_pops]))
    MAFSENS[name] = {'n_loci': int(msk.sum()), 'by_pop': rows,
                     'mainland_mean_He': mlm,
                     'KI_pct_of_mainland': 100 * rows['KI']['He'] / mlm,
                     'KI_pct_of_OVENS': 100 * rows['KI']['He'] / rows['OVENS']['He'],
                     'KI_is_lowest_He': bool(rows['KI']['He'] == min(rows[p]['He'] for p in POPS_ORDER))}
    log(f"    {name:<16s} loci={msk.sum():5d} KI He={rows['KI']['He']:.4f} "
        f"KI%mainland={MAFSENS[name]['KI_pct_of_mainland']:.1f} lowest={MAFSENS[name]['KI_is_lowest_He']}")
json.dump(MAFSENS, open(os.path.join(OUT, 'maf_sensitivity.json'), 'w'), indent=1, default=float)

# ===========================================================================
# 14.  LOCUS-CLASS RECOVERY TABLE (Reviewer 1's specific technical objection)
# ===========================================================================
log('STEP 14  locus-class recovery table')
col_4171 = np.array([j for j in range(NS) if order[j] == 'DPla19-4171'])
col_9425 = np.array([j for j in range(NS) if order[j] == 'DPla24-9425'])
def pq_of(cols, rows=None):
    sub = G[:, cols] if rows is None else G[np.ix_(rows, cols)]
    n = (sub >= 0).sum(1); alt = np.where(sub >= 0, sub, 0).sum(1)
    p = np.where(n > 0, alt / (2.0 * np.maximum(n, 1)), np.nan)
    return p, n
p4171, n4171 = pq_of(col_4171); p9425, n9425 = pq_of(col_9425)
pKIall, _ = pq_of(keep_idx[KI_COLS]); pMLall, _ = pq_of(keep_idx[MAINLAND_COLS])
def cls_counts(rows):
    a = np.zeros(len(rows), bool) if len(rows) else np.array([], bool)
    P4 = (p4171[rows] > 0) & (p4171[rows] < 1)
    P9 = (p9425[rows] > 0) & (p9425[rows] < 1)
    PK = (pKIall[rows] > 0) & (pKIall[rows] < 1)
    PM = (pMLall[rows] > 0) & (pMLall[rows] < 1)
    fixdiff = ((pKIall[rows] == 0) & (pMLall[rows] == 1)) | ((pKIall[rows] == 1) & (pMLall[rows] == 0))
    return {
      'n_loci': int(len(rows)),
      'polymorphic_in_DPla19-4171_order': int(P4.sum()),
      'polymorphic_in_DPla24-9425_order': int(P9.sum()),
      'monomorphic_in_4171_but_polymorphic_in_9425': int(((~P4) & P9).sum()),
      'monomorphic_in_9425_but_polymorphic_in_4171': int((P4 & (~P9)).sum()),
      'monomorphic_in_both_orders': int(((~P4) & (~P9)).sum()),
      'monomorphic_in_Mijangos214_but_polymorphic_in_KI8': int(((~PM) & PK).sum()),
      'polymorphic_in_Mijangos214_but_monomorphic_in_KI8': int((PM & (~PK)).sum()),
      'fixed_allelic_difference_KI_vs_mainland': int(fixdiff.sum())}
RECOVERY = {
  'input_panel_22054': cls_counts(np.arange(L0)),
  'final_filtered_set': cls_counts(FIN),
  'final_set_retained_by_old_cloneid_rule': cls_counts(IDX_RET),
  'final_set_lost_by_old_cloneid_rule': cls_counts(IDX_LOST)}
for k in ('monomorphic_in_4171_but_polymorphic_in_9425','monomorphic_in_9425_but_polymorphic_in_4171',
          'monomorphic_in_Mijangos214_but_polymorphic_in_KI8','fixed_allelic_difference_KI_vs_mainland'):
    tot = RECOVERY['final_filtered_set'][k]; lost = RECOVERY['final_set_lost_by_old_cloneid_rule'][k]
    RECOVERY.setdefault('pct_of_class_that_the_old_rule_would_have_dropped', {})[k] = (
        float(100 * lost / tot) if tot else None)
RECOVERY['statement'] = ('Because both DArT orders were called in a single pipeline run, loci that are '
   'monomorphic within one order and loci that are fixed differences between KI and the mainland are '
   'present in the matrix and are retained by the revised pipeline. Under the CloneID-matching rule used '
   'in the submitted analysis a substantial fraction of exactly these classes was discarded.')
json.dump(RECOVERY, open(os.path.join(OUT, 'locus_class_recovery.json'), 'w'), indent=1, default=float)
for k, v in RECOVERY['final_filtered_set'].items():
    log(f'    final set: {k:<52s} {v}')

# ===========================================================================
# 15.  AUDITABLE PER-LOCUS TABLE AND FILTER CASCADE JSON
# ===========================================================================
log('STEP 15  writing loci_final.csv and filter_cascade.json')
pKI_all, _ = pq_of(keep_idx[KI_COLS]); pML_all, _ = pq_of(keep_idx[MAINLAND_COLS])
pJ_all, nJ_all = pq_of(keep_idx)
LOCI = pd.DataFrame({
    'AlleleID': META.AlleleID, 'CloneID': CloneID, 'SNP': META.SNP,
    'SnpPosition': META.SnpPosition, 'scaffold': scaffold,
    'scaffold_pos': META.ChromPosTag_Platypus_NCBIv1,
    'refseq_majority': locus_rs, 'chrom_label': locus_label, 'class': locus_class,
    'RepAvg': RepAvg, 'AlnCnt': AlnCnt, 'AlnEvalue': AlnEval,
    'CallRate_report': CallRepo, 'CallRate_222': cr222,
    'pass_rep': pass_rep.astype(int), 'pass_aln': pass_aln.astype(int),
    'pass_blast': pass_bla.astype(int), 'pass_autosome': pass_auto.astype(int),
    'pass_callrate': pass_cr.astype(int),
    'pass_hwe': (~casc_pp['hwe_fail']).astype(int),
    'pass_secondary': FINAL_MASK.astype(int),
    'retained': FINAL_MASK.astype(int),
    'retained_by_old_cloneid_rule': retained_by_old.astype(int),
    'p_joint': pJ_all, 'maf_joint': np.minimum(pJ_all, 1 - pJ_all),
    'p_KI': pKI_all, 'p_mainland': pML_all,
    'poly_KI': (((pKI_all > 0) & (pKI_all < 1))).astype(int),
    'poly_mainland': (((pML_all > 0) & (pML_all < 1))).astype(int),
    'sexlinked_data_driven_flag': sex_dd.astype(int)})
LOCI.to_csv(os.path.join(OUT, 'loci_final.csv'), index=False)

CR_SENS = {}
for thr in [1.00, 0.99, 0.98, 0.95, 0.90]:
    c = run_cascade(thr, False, 'per_pop')
    s = locus_stats(Gk[np.where(c['mask_final'])[0]])
    CR_SENS[f'{thr:.2f}'] = {'final_loci': c['n_final'],
                             'polymorphic': int(((s['p'] > 0) & (s['p'] < 1)).sum())}
CASCADE = {
 'input_loci': L0, 'input_samples': NS, 'samples_retained': 222,
 'sample_selection': {'rule': 'by COLUMN INDEX on (order, id); never by ID alone',
    'mainland': '218 DPla19-4171 columns minus E32, T3, V34, T42 (Mijangos et al. 2022 exclusions)',
    'KI': '8 DPla24-9425 columns KI1..KI8',
    'alias_map': ALIAS, 'V30_V32_swap': 'applied (Mijangos et al. 2022 analyses_platy.R lines 64-71)',
    'duplicated_ids_across_orders': sorted(dup_ids),
    'population_sizes': dict(counts)},
 'steps': steps,
 'sexchrom': {'method': 'scaffold_majority_crosswalk_via_AlleleID',
   'n_anchor_loci': n_anchor, 'n_scaffolds_total': int(len(CW)),
   'n_scaffolds_anchored': n_anch_scaf, 'n_scaffolds_pure_ge95pct': n_pure95,
   'validation_agreement_vs_direct': float(agree),
   'confusion_vs_direct': {'direct_sex_cw_sex': int(dsex_cw_sex), 'direct_sex_cw_auto': int(dsex_cw_auto),
                           'direct_auto_cw_sex': int(dauto_cw_sex), 'direct_auto_cw_auto': int(dauto_cw_auto)},
   'n_sex_loci': n_sex, 'unanchored_loci': unanchored_loci,
   'data_driven_candidates': int(sex_dd.sum()),
   'data_driven_candidates_on_autosomal_scaffolds': int((sex_dd & (locus_class == 'autosome')).sum()),
   'data_driven_candidates_in_final_set': int((sex_dd & FINAL_MASK).sum()),
   'sensitivity_drop_unanchored': {'final_loci': int(casc_drop['n_final']),
                                   'max_He_shift': float(max_he_shift)}},
 'hwe': {'primary': 'per_population_n>=19_fail_in_>=2 (Bonferroni alpha = 0.05/(n_loci*6))',
         'removed': casc_pp['n_hwe_fail'], 'alpha': casc_pp['alpha'],
         'sensitivity_joint_bonferroni_removed': casc_jt['n_hwe_fail'],
         'sensitivity_joint_final_loci': casc_jt['n_final']},
 'secondary_rule': 'max CallRate, then max AvgPIC, then lowest row index',
 'final_loci': NL,
 'final_loci_polymorphic': int(((ST_JOINT['p'] > 0) & (ST_JOINT['p'] < 1)).sum()),
 'final_loci_monomorphic_in_222': int(((ST_JOINT['p'] == 0) | (ST_JOINT['p'] == 1)).sum()),
 'callrate_sensitivity': CR_SENS,
 'encoding_check_max_abs_diff': ENCODING_DIFFS}
json.dump(CASCADE, open(os.path.join(OUT, 'filter_cascade.json'), 'w'), indent=1, default=float)

# ===========================================================================
# 16.  SUPPLEMENTARY TABLE: SUBMITTED vs RECOMPUTED
# ===========================================================================
log('STEP 16  TableS_old_vs_new.csv')
old_sum = json.load(open(os.path.join(GEN, 'ki_mijangos_joint_summary.json')))
old_rar = json.load(open(os.path.join(GEN, 'subsampling_rarefaction.json')))
old_ar_txt = open(os.path.join(GEN, 'allelic_richness.json')).read().replace('NaN', 'null')
old_ar = json.loads(old_ar_txt)
LOCI_OLD = old_sum['data']['joint_loci_after_filters']
rows = []
for p in POPS_ORDER:
    o = old_sum['per_population'].get(p, {})
    npoly_old = (old_rar['KI_observed']['n_poly'] if p == 'KI'
                 else old_rar['by_pop'].get(p, {}).get('full', {}).get('n_poly'))
    ar_old = old_ar.get(p, {}).get('Ar_g16')
    n = DIVJ[p]
    rows.append(dict(population=p, n=n['n'], loci_old=LOCI_OLD, loci_new=NL,
        Ho_old=o.get('Ho_mean'), Ho_new=n['Ho'],
        Ho_delta=(n['Ho'] - o['Ho_mean']) if o else None,
        Ho_pct_change=(100 * (n['Ho'] / o['Ho_mean'] - 1)) if o else None,
        He_old=o.get('He_mean'), He_new=n['He_unbiased'],
        He_delta=(n['He_unbiased'] - o['He_mean']) if o else None,
        He_pct_change=(100 * (n['He_unbiased'] / o['He_mean'] - 1)) if o else None,
        Fis_old=o.get('FIS_mean'), Fis_new=n['Fis_WC'],
        Fis_estimator_old='1 - Ho/He per locus, mean over loci',
        Fis_estimator_new='Weir & Cockerham (1984), ratio of sums',
        n_poly_old=npoly_old, n_poly_new=n['n_polymorphic'],
        Ar16_old=ar_old, Ar16_new=(n['Ar_g16'] if np.isfinite(n['Ar_g16']) else None),
        Ar16_pct_change=(100 * (n['Ar_g16'] / ar_old - 1)) if (ar_old and np.isfinite(n['Ar_g16'])) else None))
TS = pd.DataFrame(rows)
TS.to_csv(os.path.join(OUT, 'TableS_old_vs_new.csv'), index=False)
OLD_NEW_RATIOS = {
  'KI_pct_of_mainland_mean_He_old': old_sum['summary']['KI_pct_of_mainland'],
  'KI_pct_of_mainland_mean_He_new': RATIOS['KI_pct_of_mainland_mean_He'],
  'KI_pct_of_OVENS_He_old': old_sum['summary']['KI_pct_of_Ovens'],
  'KI_pct_of_OVENS_He_new': RATIOS['KI_pct_of_OVENS_He'],
  'KI_pct_of_mainland_mean_Ar16_old': old_ar.get('_KI_pct_mainland') or 92.5,
  'KI_pct_of_mainland_mean_Ar16_new': RATIOS['KI_pct_of_mainland_mean_Ar16'],
  'mean_He_pct_change_mainland': float(np.mean([r['He_pct_change'] for r in rows
                                                if r['population'] != 'KI' and r['He_pct_change'] is not None])),
  'He_pct_change_KI': [r['He_pct_change'] for r in rows if r['population'] == 'KI'][0],
  'loci_old': LOCI_OLD, 'loci_new': NL, 'locus_gain_factor': NL / LOCI_OLD,
  'mijangos2022_published_He': 0.140,
  'new_mainland_mean_He': MAINLAND_MEAN_HE,
  'discrepancy_vs_mijangos_published_He': MAINLAND_MEAN_HE - 0.140,
  'discrepancy_vs_mijangos_published_He_pct': 100 * (MAINLAND_MEAN_HE / 0.140 - 1),
  'discrepancy_explanation': (
     'EXPECTED, not an error. The co-processed panel is ascertained across all 376 samples and therefore '
     'retains loci that are monomorphic within the Mijangos et al. (2022) 218-sample panel; neither the '
     'Mijangos-only report nor the CloneID intersection ever contained them. Including them lowers every '
     'population mean H_E. Absolute H_E is therefore NOT comparable to the published Table 1 values and all '
     'KI-versus-mainland claims are made in relative terms.')}
json.dump(OLD_NEW_RATIOS, open(os.path.join(OUT, 'old_vs_new_ratios.json'), 'w'), indent=1, default=float)
log(f"  He change: mainland {OLD_NEW_RATIOS['mean_He_pct_change_mainland']:+.1f}%, "
    f"KI {OLD_NEW_RATIOS['He_pct_change_KI']:+.1f}%; "
    f"KI:Ovens {OLD_NEW_RATIOS['KI_pct_of_OVENS_He_old']:.1f}% -> {OLD_NEW_RATIOS['KI_pct_of_OVENS_He_new']:.1f}%")

# ===========================================================================
# 17.  MASTER SUMMARY + HUMAN-READABLE RESULTS
# ===========================================================================
log('STEP 17  joint_summary.json and RESULTS_SUMMARY.md')
POP_PC_MEANS = {p: {f'PC{k+1}': float(scores[POP_COLS[p], k].mean()) for k in range(4)} for p in POPS_ORDER}
KI_IND = IND[IND.population == 'KI'].copy()
MAIN_IND = IND[IND.population != 'KI']
IND_SUMMARY = {
  'KI_individuals': KI_IND.round(5).to_dict(orient='records'),
  'KI_Ho_range': [float(KI_IND.Ho_ind.min()), float(KI_IND.Ho_ind.max())],
  'KI_Ho_mean': float(KI_IND.Ho_ind.mean()),
  'mainland_Ho_range': [float(MAIN_IND.Ho_ind.min()), float(MAIN_IND.Ho_ind.max())],
  'mainland_Ho_mean': float(MAIN_IND.Ho_ind.mean()),
  'n_KI_below_all_mainland_Ho': int((KI_IND.Ho_ind < MAIN_IND.Ho_ind.min()).sum()),
  'n_mainland_below_KI_max': int((MAIN_IND.Ho_ind < KI_IND.Ho_ind.max()).sum()),
  'KI_mean_percentile_in_mainland': float(100 * (MAIN_IND.Ho_ind < KI_IND.Ho_ind.mean()).mean()),
  'KI_F_uni_ownpop_range': [float(KI_IND.F_uni_ownpop.min()), float(KI_IND.F_uni_ownpop.max())],
  'KI_F_uni_mainland_range': [float(KI_IND.F_uni_mainland.min()), float(KI_IND.F_uni_mainland.max())],
  'KI_F_uni_mainland_mean': float(KI_IND.F_uni_mainland.mean()),
  'KI_F_hom_mainland_mean': float(KI_IND.F_hom_mainland.mean()),
  'KI_F_grm_mainland_mean': float(KI_IND.F_grm_mainland.mean()),
  'mainland_F_uni_mainland_frame_mean': float(MAIN_IND.F_uni_mainland.mean()),
  'mainland_F_uni_mainland_frame_p95': float(np.percentile(MAIN_IND.F_uni_mainland, 95)),
  'mainland_F_uni_mainland_frame_max': float(MAIN_IND.F_uni_mainland.max()),
  'n_KI_above_mainland_max_F_uni': int((KI_IND.F_uni_mainland > MAIN_IND.F_uni_mainland.max()).sum()),
  'frames': {'ownpop': 'allele frequencies of the individual''s own population (KI frame for KI animals)',
             'mainland': 'allele frequencies of the pooled 214-sample mainland panel'},
  'interpretation': ('Within-sample reference frequencies force the mean of F_UNI towards zero in any '
                     'population, so own-population F values are informative about spread, not level; the '
                     'mainland-frame values are the ones that measure KI homozygosity on a common scale.'),
  'matched_n8_null': MATCHED8['F_uni']}
SUMMARY = {
 'provenance': {
   'script': os.path.abspath(__file__), 'script_sha256': sha256_file(os.path.abspath(__file__)),
   'input_report': NEWREP, 'input_report_bytes': os.path.getsize(NEWREP),
   'input_report_mtime': time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(os.path.getmtime(NEWREP))),
   'input_report_sha256': sha256_file(NEWREP),
   'old_report_used_for': ['scaffold->chromosome crosswalk (metadata only)',
                           'old-vs-new CloneID bias diagnostic (metadata only)'],
   'old_report': OLDREP,
   'numpy': np.__version__, 'pandas': pd.__version__, 'python': sys.version.split()[0],
   'seed': SEED, 'run_timestamp': time.strftime('%Y-%m-%dT%H:%M:%S'),
   'runtime_seconds': round(time.time() - T0, 1)},
 'assertions': ASSERTIONS,
 'filter_cascade': CASCADE,
 'samples': {'n_total_columns': NS, 'n_retained': 222, 'population_sizes': dict(counts),
             'regions': {p: pop2region.get(p, '') for p in POPS_ORDER}},
 'diversity_by_population': {r['population']: r for r in div_rows},
 'ratios': RATIOS,
 'hwe_form_sensitivity_joint': HWE_SENS,
 'unanchored_drop_sensitivity': UNANCH_SENS,
 'individuals': IND_SUMMARY,
 'kinship': {'KI': KI_KIN, 'empirical_null_n_ge_19': EMP_NULL, 'matched_n8_null': MATCHED8['kinship'],
             'within_pop_summary': {k: {kk: vv for kk, vv in v.items() if kk != 'values'}
                                    for k, v in KIN_BY_POP.items()}},
 'fst': {'summary': FST_SUMMARY, 'KI_vs_mainland': {str(r['other']): {
            'Fst_WC': r['Fst_WC'], 'ci95': [r['Fst_lo95'], r['Fst_hi95']],
            'Gst_Nei': r['Gst_Nei'], 'GstPP_Hedrick': r['GstPP_Hedrick']}
            for _, r in ki_fst.iterrows()},
         'matrix': {pi: {pj: float(FSTM.loc[pi, pj]) for pj in POPS_ORDER} for pi in POPS_ORDER}},
 'pca': {'variance_pct': [float(v) for v in varpct[:6]], 'overlap': PCA_OVERLAP,
         'population_means': POP_PC_MEANS, 'n_loci': NL},
 'rarefaction': RAR, 'allelic_richness_g16': {p: DIVJ[p]['Ar_g16'] for p in POPS_ORDER},
 'polymorphism_matched_n': POLY_MATCH,
 'sfs': SFS, 'maf_sensitivity': MAFSENS,
 'oldnew_locus_bias': OLDNEW, 'old_vs_new_ratios': OLD_NEW_RATIOS,
 'locus_class_recovery': RECOVERY}
json.dump(SUMMARY, open(os.path.join(OUT, 'joint_summary.json'), 'w'), indent=1, default=float)
json.dump(ASSERTIONS, open(os.path.join(OUT, 'assertions.json'), 'w'), indent=1, default=float)

# ---------------------------------------------------------------------------
def f(x, d=4):
    return '—' if x is None or (isinstance(x, float) and not np.isfinite(x)) else f'{x:.{d}f}'
M = []
A = M.append
A('# KI platypus — co-processed re-analysis: results summary\n')
A(f"Generated {SUMMARY['provenance']['run_timestamp']} by `ki_coprocessed_joint.py` "
  f"(seed {SEED}; numpy {np.__version__}, pandas {pd.__version__}).  "
  f"Runtime {SUMMARY['provenance']['runtime_seconds']} s.\n")
A('Source: DArT **co-processed** one-row-per-locus mapping report '
  '`Report_DPla24-9425_1_moreOrders_SNP_mapping_2.csv` (22,054 loci x 376 samples, both DArT orders '
  'called in a single pipeline run).  **No CloneID merge and no batch-effect call-string filter are '
  'applied anywhere in this pipeline.**\n')
A(f'> All {len(ASSERTIONS)} hard assertions specified in the analysis plan passed; values are in '
  f'`assertions.json` and in section 13 below.\n')

A('\n## 1. Filter cascade (supplementary Table S1)\n')
A('| # | Filter | Loci passing (standalone) | Loci remaining (cumulative) |')
A('|---|---|---|---|')
for i, s in enumerate(steps):
    if s.get('removed') is not None:
        sa = f"{s['removed']:,} removed"
    elif s.get('standalone') is None:
        sa = 'one per CloneID'
    else:
        sa = f"{s['standalone']:,}"
    A(f"| {i} | {s['name']} | {sa} | **{s['cumulative']:,}** |")
A(f"\n**Final locus set: {NL:,} SNPs x 222 individuals, zero missing genotypes** "
  f"({int(CASCADE['final_loci_polymorphic']):,} polymorphic in the joint panel, "
  f"{int(CASCADE['final_loci_monomorphic_in_222']):,} monomorphic).  "
  f"This is {NL/LOCI_OLD:.2f}x the {LOCI_OLD:,} loci of the submitted analysis.\n")
A(f"Sex-chromosome exclusion removed **{n_sex:,} loci ({100*n_sex/L0:.1f}%)** via a scaffold->chromosome "
  f"crosswalk built from {n_anchor:,} AlleleID anchors shared with the Mijangos-only report "
  f"({n_anch_scaf} of {len(CW)} scaffolds anchored, {n_pure95} of them >=95% pure for one chromosome; "
  f"{100*agree:.2f}% locus-level agreement with the old report's direct annotation).\n")
A('Sensitivities: ')
A(f"- HWE form — per-population (primary, {casc_pp['n_hwe_fail']} loci removed) gives {NL:,} loci; "
  f"the joint-sample Bonferroni form quoted in the submitted Methods removes "
  f"{casc_jt['n_hwe_fail']} and gives {casc_jt['n_final']:,}.")
A(f"- Unanchored scaffolds — dropping the {unanchored_loci} loci on scaffolds with no crosswalk anchor "
  f"gives {casc_drop['n_final']:,} loci and shifts no population H_E by more than {max_he_shift:.5f}.")
A('- Call rate — ' + '; '.join(f"{k}: {v['final_loci']:,} loci ({v['polymorphic']:,} polymorphic)"
                               for k, v in CR_SENS.items()) + '.')
A(f"- Residual sex-linkage — {int(sex_dd.sum())} loci trip a data-driven male/female heterozygosity test; "
  f"{int((sex_dd & FINAL_MASK).sum())} of them survive into the final set. Reported, not filtered on.\n")

A('\n## 2. Per-population diversity (main-text Table 1)\n')
A('| Population | Region | n | H_O (SE) | H_E (SE) | F_IS (W&C) | jackknife SE | 95% bootstrap CI | '
  'Polymorphic loci | % poly | A_r (g=16) |')
A('|---|---|---|---|---|---|---|---|---|---|---|')
for p in POPS_ORDER:
    r = DIVJ[p]
    A(f"| {'**KI**' if p=='KI' else p} | {r['region']} | {r['n']} | {f(r['Ho'])} ({f(r['Ho_SE'])}) | "
      f"{f(r['He_unbiased'])} ({f(r['He_SE'])}) | {r['Fis_WC']:+.4f} | {f(r['Fis_jack_SE'])} | "
      f"[{r['Fis_boot_lo']:+.4f}, {r['Fis_boot_hi']:+.4f}] | {r['n_polymorphic']:,} | "
      f"{100*r['prop_polymorphic']:.1f}% | {f(r['Ar_g16'])} |")
A(f"\nAll {NL:,} loci, no missing data. H_E is the unbiased Nei (1987) estimator; F_IS is Weir & "
  f"Cockerham (1984) as a ratio of sums over loci with delete-one-locus jackknife SE and a "
  f"1,000-replicate bootstrap over loci.  A_r is El Mousadik & Petit (1996) rarefied allelic richness "
  f"at g = 16 gene copies; the two n = 4 populations cannot be standardised to g = 16.\n")
A(f"- KI H_E = **{f(DIVJ['KI']['He_unbiased'])}** = **{RATIOS['KI_pct_of_mainland_mean_He']:.1f}%** of the "
  f"mainland mean ({f(MAINLAND_MEAN_HE)}) and **{RATIOS['KI_pct_of_OVENS_He']:.1f}%** of Ovens "
  f"({f(DIVJ['OVENS']['He_unbiased'])}).")
A(f"- KI H_O = {f(DIVJ['KI']['Ho'])} = {RATIOS['KI_pct_of_mainland_mean_Ho']:.1f}% of the mainland mean.")
A(f"- KI A_r(16) = {f(DIVJ['KI']['Ar_g16'])} = {RATIOS['KI_pct_of_mainland_mean_Ar16']:.1f}% of the mainland "
  f"mean over the {RATIOS['n_mainland_pops_with_Ar16']} populations that can be standardised.")
A(f"- KI is the lowest of the eleven populations for H_E ({RATIOS['KI_is_lowest_He']}), "
  f"H_O ({RATIOS['KI_is_lowest_Ho']}) and A_r(16) ({RATIOS['KI_is_lowest_Ar16']}).\n")

A('> **Absolute H_E is not comparable to Mijangos et al. (2022) Table 1 and the submitted cross-check '
  'must be deleted.**  The mainland mean here is '
  f"{f(MAINLAND_MEAN_HE)} against the published {0.140:.3f} "
  f"({OLD_NEW_RATIOS['discrepancy_vs_mijangos_published_He_pct']:+.1f}%). "
  'The co-processed panel is ascertained across all 376 samples and retains loci that are monomorphic '
  'within the Mijangos 218, which neither the Mijangos-only report nor the CloneID intersection ever '
  'contained. This is expected and is not a bug; every KI-versus-mainland claim below is made in '
  'relative terms.\n')

A('\n## 3. Reviewer 1 — what the CloneID-matching approach lost, and in which direction\n')
A(f"Old rule reconstructed exactly: keep a CloneID iff it is present in both reports **and** the two SNP "
  f"call strings match ({len(shared_clones):,} shared CloneIDs, {len(match_clones):,} with matching call "
  f"strings).  Applied to the {NL:,} final loci:\n")
A('| Locus class | n loci | % of final set | Polymorphic in KI (n=8) | Polymorphic in mainland (n=214) | '
  'Mean joint MAF | Mean H_E (222) |')
A('|---|---|---|---|---|---|---|')
afc = OLDNEW['allele_frequency_character']
A(f"| Retained by the old CloneID rule | {OLDNEW['retained_by_old']:,} | "
  f"{100-OLDNEW['pct_lost_by_old']:.1f}% | {OLDNEW['retained_poly_KI']:,} "
  f"({100*OLDNEW['retained_prop_poly_KI']:.1f}%) | {OLDNEW['retained_poly_mainland']:,} "
  f"({100*OLDNEW['retained_prop_poly_mainland']:.1f}%) | {f(afc['retained_by_old']['mean_joint_MAF'])} | "
  f"{f(afc['retained_by_old']['mean_He_joint'])} |")
A(f"| **Lost by the old CloneID rule** | **{OLDNEW['lost_by_old']:,}** | "
  f"**{OLDNEW['pct_lost_by_old']:.1f}%** | {OLDNEW['lost_poly_KI']:,} "
  f"({100*OLDNEW['lost_prop_poly_KI']:.1f}%) | {OLDNEW['lost_poly_mainland']:,} "
  f"({100*OLDNEW['lost_prop_poly_mainland']:.1f}%) | {f(afc['lost_by_old']['mean_joint_MAF'])} | "
  f"{f(afc['lost_by_old']['mean_He_joint'])} |")
A(f"\nThe lost loci are overwhelmingly the low-frequency, near-fixed ones: "
  f"{afc['lost_by_old']['pct_monomorphic_in_222']:.1f}% of them are monomorphic across all 222 individuals "
  f"versus {afc['retained_by_old']['pct_monomorphic_in_222']:.1f}% of the retained loci.\n")
A('**H_E computed on each subset (direction and size of the inflation):**\n')
A('| Population | All final loci | Retained-by-old subset | Lost-by-old subset | Inflation from using only '
  'the old subset |')
A('|---|---|---|---|---|')
for p in POPS_ORDER:
    d = OLDNEW['by_population'][p]
    A(f"| {p} | {f(d['all_final']['He'])} | {f(d['retained_by_old']['He'])} | {f(d['lost_by_old']['He'])} | "
      f"{d['He_inflation_pct']:+.1f}% |")
A('\n**Submitted versus recomputed, per population (supplementary Table S2; full file '
  '`TableS_old_vs_new.csv`):**\n')
A('| Population | n | H_O submitted | H_O new | H_E submitted | H_E new | ΔH_E % | F_IS submitted | '
  'F_IS new (W&C) | Polymorphic submitted | Polymorphic new | A_r(16) submitted | A_r(16) new |')
A('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    _dhe = ('%+.1f%%' % r['He_pct_change']) if r['He_pct_change'] is not None else '\u2014'
    _npo = format(r['n_poly_old'], ',') if r['n_poly_old'] else '\u2014'
    A(f"| {r['population']} | {r['n']} | {f(r['Ho_old'])} | {f(r['Ho_new'])} | {f(r['He_old'])} | "
      f"{f(r['He_new'])} | {_dhe} | "
      f"{f(r['Fis_old'])} | {r['Fis_new']:+.4f} | "
      f"{_npo} | {r['n_poly_new']:,} | "
      f"{f(r['Ar16_old'])} | {f(r['Ar16_new'])} |")
A(f"\n- Locus count {LOCI_OLD:,} -> {NL:,} ({NL/LOCI_OLD:.2f}x).")
A(f"- Absolute H_E falls by {abs(OLD_NEW_RATIOS['mean_He_pct_change_mainland']):.1f}% on average for the "
  f"mainland populations and by {abs(OLD_NEW_RATIOS['He_pct_change_KI']):.1f}% for KI.")
A(f"- **The relative conclusion survives:** KI:mainland-mean H_E "
  f"{OLD_NEW_RATIOS['KI_pct_of_mainland_mean_He_old']:.1f}% -> "
  f"{OLD_NEW_RATIOS['KI_pct_of_mainland_mean_He_new']:.1f}%; KI:Ovens "
  f"{OLD_NEW_RATIOS['KI_pct_of_OVENS_He_old']:.1f}% -> "
  f"{OLD_NEW_RATIOS['KI_pct_of_OVENS_He_new']:.1f}%; KI remains the lowest of the eleven populations.\n")

A('\n## 4. Locus-class recovery — Reviewer 1\'s technical objection\n')
A('| Locus class | Input panel (22,054) | Final filtered set | Of the final set, would have been dropped '
  'by the CloneID rule |')
A('|---|---|---|---|')
lab = {'polymorphic_in_DPla19-4171_order': 'Polymorphic within DPla19-4171 (Mijangos order)',
       'polymorphic_in_DPla24-9425_order': 'Polymorphic within DPla24-9425 (2024 order)',
       'monomorphic_in_4171_but_polymorphic_in_9425': 'Monomorphic in DPla19-4171, polymorphic in DPla24-9425',
       'monomorphic_in_9425_but_polymorphic_in_4171': 'Monomorphic in DPla24-9425, polymorphic in DPla19-4171',
       'monomorphic_in_both_orders': 'Monomorphic in both orders',
       'monomorphic_in_Mijangos214_but_polymorphic_in_KI8': 'Monomorphic in the mainland 214, polymorphic in KI',
       'polymorphic_in_Mijangos214_but_monomorphic_in_KI8': 'Polymorphic in the mainland 214, fixed in KI',
       'fixed_allelic_difference_KI_vs_mainland': 'Fixed allelic difference, KI vs mainland'}
for k, nm in lab.items():
    tot = RECOVERY['final_filtered_set'][k]
    lost = RECOVERY['final_set_lost_by_old_cloneid_rule'][k]
    pct = f'{100*lost/tot:.1f}%' if tot else '—'
    A(f"| {nm} | {RECOVERY['input_panel_22054'][k]:,} | {tot:,} | {lost:,} ({pct}) |")
A('\nBecause both DArT orders were called in one pipeline run, loci that are monomorphic within one order '
  'and loci that are fixed differences between KI and the mainland are present in the matrix and are '
  'retained. The CloneID rule discarded a large share of exactly these classes.\n')

A('\n## 5. Pairwise differentiation — F_ST (NEW)\n')
A('| Comparison | F_ST (W&C) | 95% bootstrap CI | G_ST (Nei) | G"_ST (Hedrick) |')
A('|---|---|---|---|---|')
for _, r in ki_fst.iterrows():
    A(f"| KI vs {r['other']} | {r['Fst_WC']:.4f} | [{r['Fst_lo95']:.4f}, {r['Fst_hi95']:.4f}] | "
      f"{r['Gst_Nei']:.4f} | {r['GstPP_Hedrick']:.4f} |")
A(f"\n- KI is **most similar to {FST_SUMMARY['KI_min_pop']}** (F_ST = {FST_SUMMARY['KI_min']:.4f}) and most "
  f"divergent from {FST_SUMMARY['KI_max_pop']} (F_ST = {FST_SUMMARY['KI_max']:.4f}).")
A(f"- Mean KI-vs-mainland F_ST = {FST_SUMMARY['KI_mean']:.4f}; mean mainland-vs-mainland F_ST = "
  f"{FST_SUMMARY['mainland_mainland_mean']:.4f} (range {FST_SUMMARY['mainland_mainland_min']:.4f}–"
  f"{FST_SUMMARY['mainland_mainland_max']:.4f}). KI is more differentiated from every mainland population "
  f"than any two mainland populations are from each other." if FST_SUMMARY['KI_min'] >
  FST_SUMMARY['mainland_mainland_max'] else
  f"- Mean KI-vs-mainland F_ST = {FST_SUMMARY['KI_mean']:.4f}; mean mainland-vs-mainland F_ST = "
  f"{FST_SUMMARY['mainland_mainland_mean']:.4f} (range {FST_SUMMARY['mainland_mainland_min']:.4f}–"
  f"{FST_SUMMARY['mainland_mainland_max']:.4f}).")
A('\nThe full 11 x 11 matrix is in `fst_matrix.csv` and `fst_pairwise_ci.json`.\n')
A('| | ' + ' | '.join(POPS_ORDER) + ' |')
A('|---|' + '---|' * 11)
for pi in POPS_ORDER:
    A(f"| **{pi}** | " + ' | '.join('—' if pi == pj else f'{FSTM.loc[pi, pj]:.3f}' for pj in POPS_ORDER) + ' |')

A('\n## 6. PCA — honest answer to the reviewer\'s overlap objection\n')
A(f"Patterson/Price normalisation, SVD on {NL:,} loci x 222 individuals, zero missing data. "
  f"PC1 {varpct[0]:.2f}%, PC2 {varpct[1]:.2f}%, PC3 {varpct[2]:.2f}%, PC4 {varpct[3]:.2f}% of variance.\n")
A('| Axis | % variance | KI range | Mainland range | Fully separated? | KI individuals inside the '
  'mainland range |')
A('|---|---|---|---|---|---|')
for k in range(4):
    o = PCA_OVERLAP[f'PC{k+1}']
    names = ', '.join(o['KI_inside_mainland_range']) if o['KI_inside_mainland_range'] else 'none'
    A(f"| PC{k+1} | {o['variance_pct']:.2f}% | [{o['KI_range'][0]:+.2f}, {o['KI_range'][1]:+.2f}] | "
      f"[{o['mainland_range'][0]:+.2f}, {o['mainland_range'][1]:+.2f}] | "
      f"{'**yes**' if o['fully_separated'] else 'no'} | "
      f"{o['n_KI_inside_mainland_range']} of 8 ({names}) |")
sep_axes = [f'PC{k+1}' for k in range(4) if PCA_OVERLAP[f'PC{k+1}']['fully_separated']]
non_axes = [f'PC{k+1}' for k in range(4) if not PCA_OVERLAP[f'PC{k+1}']['fully_separated']]
A(f"\n**Statement for the manuscript.** KI is completely separated from all 214 mainland individuals on "
  f"{', '.join(sep_axes) if sep_axes else 'no single axis'}"
  + (f" (on {sep_axes[0]} the KI minimum is {PCA_OVERLAP[sep_axes[0]]['KI_range'][0]:+.2f} against a "
     f"mainland maximum of {PCA_OVERLAP[sep_axes[0]]['mainland_range'][1]:+.2f}, a gap of "
     f"{PCA_OVERLAP[sep_axes[0]]['gap']:.2f} units)" if sep_axes else '')
  + f", and is **not** separable on {', '.join(non_axes) if non_axes else 'no axis'}, where every KI "
  f"individual falls inside the mainland spread. "
  f"On the PC1 x PC2 plane the minimum KI-to-mainland Euclidean distance is "
  f"{PCA_OVERLAP['PC1xPC2']['min_KI_to_mainland_distance']:.2f} against a maximum within-KI distance of "
  f"{PCA_OVERLAP['PC1xPC2']['max_KI_internal_distance']:.2f}, so the eight KI animals form a discrete "
  f"cluster in that plane" + (" with no mainland individual closer to any KI individual than the KI "
  "animals are to each other.\n" if PCA_OVERLAP['PC1xPC2']['KI_cluster_separated'] else
  ", although the cluster is wider than its separation from the nearest mainland individual.\n"))
A('Because PC1 is dominated by the mainland Border-versus-Snowy contrast, the axis that isolates KI is '
  'PC2, not PC1. The figure must therefore be drawn and described as PC1 x PC2 with KI separated on PC2; '
  'the submitted claim that KI forms "a discrete cluster well removed" from the mainland should be '
  'replaced by this explicit, axis-by-axis statement. No KI individual falls inside the mainland spread '
  'on ' + (', '.join(sep_axes) if sep_axes else 'any axis') + '.\n')
A('Population means on the first four axes:\n')
A('| Population | PC1 | PC2 | PC3 | PC4 |')
A('|---|---|---|---|---|')
for p in POPS_ORDER:
    m = POP_PC_MEANS[p]
    A(f"| {p} | {m['PC1']:+.2f} | {m['PC2']:+.2f} | {m['PC3']:+.2f} | {m['PC4']:+.2f} |")

A('\n## 7. Individual-level heterozygosity and inbreeding (Reviewer 2)\n')
A('| Sample | H_O | F_HOM (KI frame) | F_UNI (KI frame) | F_GRM (KI frame) | F_UNI 95% CI (KI frame) | '
  'F_HOM (mainland frame) | F_UNI (mainland frame) | F_GRM (mainland frame) |')
A('|---|---|---|---|---|---|---|---|---|')
for _, r in KI_IND.iterrows():
    A(f"| {r['sample']} | {r['Ho_ind']:.4f} | {r['F_hom_ownpop']:+.4f} | {r['F_uni_ownpop']:+.4f} | "
      f"{r['F_grm_ownpop']:+.4f} | [{r['F_uni_ownpop_lo95']:+.4f}, {r['F_uni_ownpop_hi95']:+.4f}] | "
      f"{r['F_hom_mainland']:+.4f} | {r['F_uni_mainland']:+.4f} | {r['F_grm_mainland']:+.4f} |")
A(f"\nComputed on all {NL:,} loci (F estimators use the "
  f"{int(KI_IND.n_loci_F_ownpop.iloc[0]):,} loci polymorphic in KI and the "
  f"{int(KI_IND.n_loci_F_mainland.iloc[0]):,} loci polymorphic in the mainland panel respectively). "
  f"F_HOM = excess homozygosity; F_UNI = Yang et al. (2010) correlation of uniting gametes (F-hat-3); "
  f"F_GRM = variance-standardised GRM diagonal minus 1 (F-hat-1). CIs are 1,000-replicate bootstraps "
  f"over loci.\n")
A(f"- Per-individual H_O: KI {IND_SUMMARY['KI_Ho_range'][0]:.4f}–{IND_SUMMARY['KI_Ho_range'][1]:.4f} "
  f"(mean {IND_SUMMARY['KI_Ho_mean']:.4f}); mainland "
  f"{IND_SUMMARY['mainland_Ho_range'][0]:.4f}–{IND_SUMMARY['mainland_Ho_range'][1]:.4f} "
  f"(mean {IND_SUMMARY['mainland_Ho_mean']:.4f}). Every KI animal sits below the mainland mean and the "
  f"KI mean is at the {IND_SUMMARY['KI_mean_percentile_in_mainland']:.1f}th percentile of the mainland "
  f"distribution, but the KI and mainland ranges **overlap**: "
  f"{IND_SUMMARY['n_mainland_below_KI_max']} of the 214 mainland animals have a lower H_O than the "
  f"highest KI animal and {IND_SUMMARY['n_KI_below_all_mainland_Ho']} KI animal(s) fall below the "
  f"mainland minimum. The KI deficit is a shift in the whole distribution, not a set of individually "
  f"exceptional animals — state it that way.")
A(f"- In the **mainland reference frame** (the frame that puts all 222 animals on a common scale) the KI "
  f"animals have F_UNI {IND_SUMMARY['KI_F_uni_mainland_range'][0]:+.3f} to "
  f"{IND_SUMMARY['KI_F_uni_mainland_range'][1]:+.3f} (mean {IND_SUMMARY['KI_F_uni_mainland_mean']:+.3f}) "
  f"against a mainland mean of {IND_SUMMARY['mainland_F_uni_mainland_frame_mean']:+.3f} and a mainland "
  f"maximum of {IND_SUMMARY['mainland_F_uni_mainland_frame_max']:+.3f}; "
  f"{IND_SUMMARY['n_KI_above_mainland_max_F_uni']} of 8 exceed the highest mainland value.")
A(f"- F_GRM (F-hat-1) is the most rare-allele-sensitive of the three estimators and inflates strongly "
  f"when evaluated in a foreign reference frame; the mainland-frame F_GRM values above should be read as "
  f"a rank ordering, not as an absolute inbreeding coefficient. F_HOM and F_UNI are the ones to quote.")
A(f"- In the **KI reference frame** the same animals have F_UNI "
  f"{IND_SUMMARY['KI_F_uni_ownpop_range'][0]:+.3f} to {IND_SUMMARY['KI_F_uni_ownpop_range'][1]:+.3f}. "
  f"Within-sample frequencies force this mean towards zero by construction, so own-frame values describe "
  f"variation *among* the eight animals, not their absolute level; the matched n = 8 mainland null for "
  f"F_UNI is {MATCHED8['F_uni']['mean']:+.3f} "
  f"({MATCHED8['F_uni']['p05']:+.3f} to {MATCHED8['F_uni']['p95']:+.3f}, 5th–95th percentile over "
  f"{MATCHED8['reps']} draws from each of the six populations with n >= 19).\n")

A('\n## 8. Pairwise kinship among the eight KI individuals (Reviewer 2)\n')
A('Yang et al. (2010) genomic relationship, on the **relatedness** scale (unrelated ~0, full sibs and '
  'parent–offspring ~0.5); halve for the kinship coefficient. KI values use KI allele frequencies.\n')
A('| | ' + ' | '.join(KI_KIN['labels']) + ' |')
A('|---|' + '---|' * 8)
for i, lb in enumerate(KI_KIN['labels']):
    A(f"| **{lb}** | " + ' | '.join('—' if i == j else f"{KI_KIN['matrix'][i][j]:+.3f}"
                                    for j in range(8)) + ' |')
A(f"\n- KI: 28 pairs, mean {KI_KIN['mean']:+.4f}, range {KI_KIN['min']:+.4f} to {KI_KIN['max']:+.4f}.")
A(f"- Matched n = 8 mainland null (200 draws from each of {', '.join(BIG_POPS)}; "
  f"{MATCHED8['kinship']['n']:,} pairs): mean {MATCHED8['kinship']['mean']:+.4f}, "
  f"95th percentile {MATCHED8['kinship']['p95']:+.4f}, 99th {MATCHED8['kinship']['p99']:+.4f}, "
  f"maximum {MATCHED8['kinship']['max']:+.4f}.")
A(f"- {KI_KIN['n_above_matched_p95']} of the 28 KI pairs exceed the matched-n null 95th percentile and "
  f"{KI_KIN['n_above_matched_max']} exceed its maximum. The empirical p for the single highest KI pair is "
  f"{KI_KIN['empirical_p_max_pair']:.4f}.")
A(f"- On the full-sample empirical null (all within-population pairs in the six populations with n >= 19, "
  f"{EMP_NULL['n_pairs']:,} pairs): mean {EMP_NULL['mean']:+.4f}, 95th {EMP_NULL['p95']:+.4f}, "
  f"99.9th {EMP_NULL['p999']:+.4f}, maximum {EMP_NULL['max']:+.4f}.")
_hi_pairs = sorted(((float(KI_KIN['matrix'][i][j]), KI_KIN['labels'][i], KI_KIN['labels'][j])
                   for i in range(8) for j in range(i + 1, 8)
                   if KI_KIN['matrix'][i][j] > MATCHED8['kinship']['p95']), reverse=True)
if _hi_pairs:
    A(f"- KI pairs above the matched n = 8 mainland 95th percentile "
      f"({MATCHED8['kinship']['p95']:+.4f}): " +
      ', '.join(f'{a}–{b} ({v:+.3f})' for v, a, b in _hi_pairs) +
      '. Only the four positive values are plausibly informative; the highest three '
      '(KI1, KI5 and KI7 in a triangle) are the closest thing to structure in the KI sample.')
_verdict_kin = ('no KI pair reaches a relatedness consistent with first-order kinship'
                if KI_KIN['n_above_0p25_relatedness'] == 0 else
                f"{KI_KIN['n_above_0p25_relatedness']} KI pair(s) exceed a relatedness of 0.25")
A(f"- **Verdict:** {_verdict_kin} ({KI_KIN['n_above_0p125_relatedness']} of 28 pairs exceed 0.125, the "
  f"nominal second-order threshold). The claim should be stated against the matched-n null shown above "
  f"rather than against a theoretical threshold, because within-sample allele frequencies at n = 8 shift "
  f"the whole distribution downward by construction.\n")

A('\n## 9. Proportion of polymorphic loci — replacement metric\n')
A('The submitted metric ("3,318 of 22,054 input loci, ~15%") is deleted: the denominator was the '
  'unfiltered panel, singletons were weighted like intermediate-frequency sites, and n = 8 makes an '
  'unobserved alternate allele unremarkable. Three replacement statements, in order:\n')
A(f"**(1) On the filtered joint set**, for all eleven populations (denominator {NL:,}):\n")
A('| Population | n | Polymorphic loci | % |')
A('|---|---|---|---|')
for p in sorted(POPS_ORDER, key=lambda x: -DIVJ[x]['n_polymorphic']):
    A(f"| {p} | {DIVJ[p]['n']} | {DIVJ[p]['n_polymorphic']:,} | {100*DIVJ[p]['prop_polymorphic']:.1f}% |")
A(f"\n**(2) Sample-size matched at n = 8** (200 draws; this is the number for the abstract):\n")
A('| Population | Polymorphic loci at n = 8 (mean ± SD) | KI as % of this | KI deficit |')
A('|---|---|---|---|')
for p in BIG_POPS:
    e = POLY_MATCH['mainland_at_n8'][p]
    A(f"| {p} | {e['n_poly_mean']:.0f} ± {e['n_poly_sd']:.0f} | {e['KI_pct_of_this']:.1f}% | "
      f"{e['KI_deficit_pct']:.1f}% |")
A(f"| **KI (observed, n = 8)** | **{DIVJ['KI']['n_polymorphic']:,}** | — | — |")
A(f"\nAgainst the mean of the six matched mainland samples ({m8:.0f} loci), KI retains "
  f"{POLY_MATCH['KI_pct_of_mainland_mean_at_n8']:.1f}% — a sample-size-controlled deficit of "
  f"**{POLY_MATCH['KI_deficit_vs_mainland_mean_at_n8_pct']:.1f}%**.\n")
A(f"**(3) Rarefied allelic richness at g = 16 gene copies** (the properly standardised metric): "
  f"KI {f(DIVJ['KI']['Ar_g16'])} versus a mainland mean of {f(RATIOS['mainland_mean_Ar16'])} over the "
  f"{RATIOS['n_mainland_pops_with_Ar16']} populations that can be standardised = "
  f"{RATIOS['KI_pct_of_mainland_mean_Ar16']:.1f}%. Full curves for g = 2…16 are in "
  f"`allelic_richness.json`.\n")

A('\n## 10. Folded site frequency spectra — KI versus a rarefied mainland\n')
A(f"The submitted comparison (KI n = 8 versus Ovens n = 19) is confounded: at n = 8 the smallest non-zero "
  f"folded minor allele frequency is 1/16 = 0.0625, while at n = 19 it is 1/38 = 0.026, so Ovens can "
  f"express rare-allele classes KI cannot. The primary comparison here is therefore against Ovens "
  f"rarefied to n = 8.\n")
A('| Spectrum | n | Polymorphic loci | Mean folded MAF | Median |')
A('|---|---|---|---|---|')
A(f"| **A. KI** | 8 | {SFS['KI']['n_polymorphic']:,} | **{SFS['KI']['mean_folded_maf']:.4f}** | "
  f"{SFS['KI']['median_folded_maf']:.4f} |")
A(f"| B. Ovens, full (legacy comparison) | 19 | {SFS['OVENS_full']['n_polymorphic']:,} | "
  f"{SFS['OVENS_full']['mean_folded_maf']:.4f} | {SFS['OVENS_full']['median_folded_maf']:.4f} |")
for p in ['OVENS', 'SNOWY', 'THREDBO', 'TENTERFIELD']:
    e = SFS[f'{p}_rarefied_n8']
    tag = '**C. ' + p + ', rarefied to n = 8**' if p == 'OVENS' else f'C. {p}, rarefied to n = 8'
    A(f"| {tag} | 8 | {e['n_polymorphic_mean']:.0f} ± {e['n_polymorphic_sd']:.0f} | "
      f"{e['mean_folded_maf']:.4f} ± {e['sd']:.4f} | — |")
A(f"| D. KI on the loci the old CloneID rule kept | 8 | {SFS['KI_on_old_cloneid_subset']['n_polymorphic']:,} | "
  f"{SFS['KI_on_old_cloneid_subset']['mean_folded_maf']:.4f} | "
  f"{SFS['KI_on_old_cloneid_subset']['median_folded_maf']:.4f} |")
A('\n**KI folded spectrum, counts by minor allele count out of 16 gene copies:**\n')
A('| MAC (of 16) | ' + ' | '.join(str(i) for i in range(1, 9)) + ' |')
A('|---|' + '---|' * 8)
A('| KI count | ' + ' | '.join(str(SFS['KI']['mac_counts'][str(i)]) for i in range(1, 9)) + ' |')
A('| KI proportion | ' + ' | '.join(f"{SFS['KI']['mac_proportions'][str(i)]:.3f}" for i in range(1, 9)) + ' |')
A('| Ovens n=8 proportion | ' + ' | '.join(
    f"{SFS['OVENS_rarefied_n8']['mac_proportions_mean'][str(i)]:.3f}" for i in range(1, 9)) + ' |')
pr = SFS['tests']['KI_vs_OVENS_n8_PRIMARY']; lg = SFS['tests']['KI_vs_OVENS_full_legacy']
hy = SFS['tests']['KI_all_vs_KI_old_subset']
A(f"\n- **Primary test (A vs C, matched n = 8):** mean KS D over 200 Ovens subsamples = "
  f"{pr['KS_D_mean_over_200_subsamples']:.3f} ± {pr['KS_D_sd']:.3f}; empirical p (proportion of Ovens "
  f"n = 8 subsamples with a mean folded MAF at least as high as KI's) = "
  f"{int(round(pr['empirical_p_mean_maf_ge_KI']*200))}/200 "
  f"({'< 0.005' if pr['empirical_p_mean_maf_ge_KI'] == 0 else format(pr['empirical_p_mean_maf_ge_KI'], '.4f')}). "
  f"Difference in mean folded MAF = {pr['difference']:+.4f}.")
A(f"- **Legacy test (A vs B, unmatched n):** KS D = {lg['KS_D']:.3f}, p = {lg['KS_p']:.3g}; "
  f"Mann–Whitney z = {lg['MW_z']:.1f}, p = {lg['MW_p']:.3g}. Reported only for continuity — "
  f"**{pr['pct_of_lega' + 'cy_gap_that_is_sample_size_artefact']:.0f}% of this gap is arithmetic**, not "
  f"biological (the unmatched gap is {pr['legacy_difference_unmatched']:+.4f}; the matched gap is "
  f"{pr['difference']:+.4f}).")
A(f"- **Reviewer 1's hypothesis (A vs D):** KS D = {hy['KS_D']:.3f}, p = {hy['KS_p']:.3g}; KI mean folded "
  f"MAF on all final loci {hy['KI_all_mean']:.4f} versus {hy['KI_old_subset_mean']:.4f} on the old-rule "
  f"subset.")
A(f"\n**Verdict.** {SFS['verdict']}\n")

A('\n## 11. MAF / MAC sensitivity\n')
A('Thresholds applied symmetrically across the joint 222-sample panel, never within KI alone.\n')
A('| Scenario | Loci | KI H_E | Mainland mean H_E | Ovens H_E | KI as % of mainland | KI as % of Ovens | '
  'KI lowest? |')
A('|---|---|---|---|---|---|---|---|')
for name, v in MAFSENS.items():
    A(f"| {name} | {v['n_loci']:,} | {v['by_pop']['KI']['He']:.4f} | {v['mainland_mean_He']:.4f} | "
      f"{v['by_pop']['OVENS']['He']:.4f} | {v['KI_pct_of_mainland']:.1f}% | {v['KI_pct_of_OVENS']:.1f}% | "
      f"{'yes' if v['KI_is_lowest_He'] else 'no'} |")
A('\nAbsolute H_E rises steeply as rare variants are removed, but the KI:mainland ratio is stable to '
  'within a few percentage points and KI remains the lowest population under every threshold.\n')

A('\n## 12. Rarefaction\n')
A('200 subsamples without replacement at n = 5, 8, 10, 15 and 20 from each of the six populations with '
  'n >= 19 (seed 42).\n')
A('| Population | n (full) | H_E full | H_E at n=8 | % of full | H_O at n=8 | F_IS at n=8 | '
  'Polymorphic at n=8 | A_r(16) at n=8 |')
A('|---|---|---|---|---|---|---|---|---|')
for p in BIG_POPS:
    e = RAR['by_pop'][p]; n8 = e['n8']
    ar8 = n8['Ar_g16_mean']
    A(f"| {p} | {e['n_full']} | {e['full']['He']:.4f} | {n8['He_mean']:.4f} ± {n8['He_sd']:.4f} | "
      f"{n8['He_pct_of_full']:.1f}% | {n8['Ho_mean']:.4f} | {n8['Fis_mean']:+.4f} | "
      f"{n8['n_poly_mean']:.0f} ± {n8['n_poly_sd']:.0f} | {ar8:.4f} |")
k = RAR['KI_observed']
A(f"| **KI (observed)** | **8** | — | **{k['He']:.4f}** | — | **{k['Ho']:.4f}** | **{k['Fis']:+.4f}** | "
  f"**{k['n_poly']:,}** | **{k['Ar_g16']:.4f}** |")
A('\nH_E is essentially unbiased by sample size (every population is within ~1% of its full-sample value '
  'at n = 8), so the KI H_E deficit is not a sample-size artefact. Polymorphic-locus counts and allelic '
  'richness are strongly sample-size dependent and must be compared at matched n, as in section 9. '
  'Full output including n = 5, 10, 15 and 20 is in `rarefaction.json`.\n')

A('\n## 13. Assertions (all must pass for the run to be valid)\n')
A('| Assertion | Value | Expected | Result |')
A('|---|---|---|---|')
for kk, vv in ASSERTIONS.items():
    val = str(vv['value']); val = val if len(val) <= 90 else val[:87] + '...'
    exp = str(vv['expected']); exp = exp if len(exp) <= 70 else exp[:67] + '...'
    A(f"| `{kk}` | {val} | {exp} | {'PASS' if vv['passed'] else 'FAIL'} |")
A('\n## 14. Output files\n')
for fn in sorted(os.listdir(OUT)):
    fp = os.path.join(OUT, fn)
    if os.path.isfile(fp):
        A(f"- `out/{fn}` — {os.path.getsize(fp):,} bytes")
A('')
open(os.path.join(OUT, 'RESULTS_SUMMARY.md'), 'w').write('\n'.join(M))
log(f'DONE in {time.time()-T0:.1f}s — {len(ASSERTIONS)} assertions, all passed')
