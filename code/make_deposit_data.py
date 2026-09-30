#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_deposit_data.py
====================
Builds the public deposit data tables for

    Bino, Hawke, Baring & Gongora
    "Drifting alone: genome-wide diversity in the isolated, introduced
     Kangaroo Island platypus (Ornithorhynchus anatinus)"

from (a) the DArT co-processed one-row-per-locus mapping report and
     (b) the analysis outputs in <Revision>/work/out/ produced by
         ki_coprocessed_joint.py.

It writes <Revision>/deposit/data/ and <Revision>/deposit/results/.

Why this script exists
----------------------
ki_coprocessed_joint.py holds the filtered 4,002 x 222 dosage matrix only in
memory; out/ stores per-locus metadata (loci_final.csv) and per-sample metadata
(samples_used.csv) but not the genotypes themselves.  This script reconstructs
the genotype matrix using exactly the locus mask (`retained == 1`) and the
sample mask (`included == 1`) that the pipeline wrote, then asserts the
reconstruction against the pipeline's own published counts.

The full 22,054-locus x 376-sample co-processed report is NOT copied anywhere
into the deposit: 218 of its samples belong to Mijangos et al. (2022) and are
not ours to redistribute.  Only the 4,002 retained loci x 222 analysed
individuals are written out.

Dependencies: python3, numpy, pandas.  Seed-free (no random numbers).
"""
import os, csv, json, shutil, sys
import numpy as np
import pandas as pd

HOME = os.path.expanduser('~')
GEN  = os.path.join(HOME, 'mnt', 'Kangaroo Island', 'Genetics')
REV  = os.path.join(GEN, 'Revision')
WORK = os.path.join(REV, 'work')
OUT  = os.path.join(WORK, 'out')
DEP  = os.path.join(REV, 'deposit')
DATA = os.path.join(DEP, 'data')
RES  = os.path.join(DEP, 'results')
NEWREP = os.path.join(HOME, 'mnt', 'Platypus - Genomics - General', 'DaRT Data',
                      'DPla24-9425', 'Report-DPla24-9425',
                      'Report_DPla24-9425_1_moreOrders_SNP_mapping_2.csv')

for d in (DATA, RES, os.path.join(DEP, 'code'), os.path.join(DEP, 'figures')):
    os.makedirs(d, exist_ok=True)

N_META = 26
DOSE   = {'0': 0, '1': 2, '2': 1}      # ref-hom -> 0, snp-hom -> 2, het -> 1


def log(m):
    print(m, flush=True)


# ---------------------------------------------------------------- masks
log('reading out/loci_final.csv and out/samples_used.csv')
LOCI = pd.read_csv(os.path.join(OUT, 'loci_final.csv'))
SAMP = pd.read_csv(os.path.join(OUT, 'samples_used.csv'))

assert len(LOCI) == 22054, len(LOCI)
keep_locus = LOCI['retained'].values.astype(bool)
assert keep_locus.sum() == 4002, keep_locus.sum()

inc = SAMP[SAMP['included'] == 1].copy()
assert len(inc) == 222, len(inc)
keep_cols = inc['col_index'].values.astype(int)
samp_names = inc['id_canonical'].astype(str).tolist()
assert len(set(samp_names)) == 222, 'duplicate sample names among the 222'

# ---------------------------------------------------------------- report
log('parsing the co-processed DArT report (22,054 x 376)')
with open(NEWREP, newline='', encoding='utf-8', errors='replace') as fh:
    rdr = csv.reader(fh)
    head = [next(rdr) for _ in range(7)]
    header = head[6]
    data_rows = [r for r in rdr if len(r) > 1]

assert len(header) == 402, len(header)
assert len(data_rows) == 22054, len(data_rows)
report_alleleids = [r[0] for r in data_rows]
assert report_alleleids == LOCI['AlleleID'].astype(str).tolist(), \
    'loci_final.csv row order does not match the report row order'

sel_rows = np.where(keep_locus)[0]
G = np.empty((len(sel_rows), 222), dtype=np.int8)
for k, i in enumerate(sel_rows):
    gt = data_rows[i][N_META:]
    G[k] = [DOSE.get(gt[j], -1) for j in keep_cols]
del data_rows

log(f'  reconstructed dosage matrix {G.shape}, missing cells = {int((G < 0).sum())}')
assert G.shape == (4002, 222), G.shape
assert int((G < 0).sum()) == 0, 'missing genotypes in a matrix the pipeline asserts is complete'

# ---------------------------------------------------------------- cross-checks
# per-individual H_O must reproduce out/diversity_by_individual.csv
IND = pd.read_csv(os.path.join(OUT, 'diversity_by_individual.csv')).set_index('sample')
ho = (G == 1).mean(axis=0)
ref = IND.loc[samp_names, 'Ho_ind'].values
d = float(np.max(np.abs(ho - ref)))
log(f'  max |H_O(reconstructed) - H_O(pipeline)| over 222 individuals = {d:.3e}')
assert d < 1e-9, d

# per-population polymorphic-locus counts must reproduce diversity_by_population.csv
POP = pd.read_csv(os.path.join(OUT, 'diversity_by_population.csv'))
popv = inc['population'].astype(str).values
for _, row in POP.iterrows():
    m = popv == row['population']
    sub = G[:, m]
    p = sub.sum(axis=0) / (2.0 * m.sum()) if False else sub.sum(axis=1) / (2.0 * m.sum())
    npoly = int(((p > 0) & (p < 1)).sum())
    assert npoly == int(row['n_polymorphic']), (row['population'], npoly, row['n_polymorphic'])
log('  per-population polymorphic-locus counts reproduce the pipeline exactly')

ki_cols = [s for s in samp_names if s.startswith('KI')]
assert ki_cols == [f'KI{i}' for i in range(1, 9)] or sorted(ki_cols) == sorted(
    [f'KI{i}' for i in range(1, 9)]), ki_cols
log(f'  KI columns present: {sorted(ki_cols, key=lambda s: int(s[2:]))}')

# ---------------------------------------------------------------- data/ tables
log('writing data/genotypes_filtered_4002x222.csv')
GT = pd.DataFrame(G, columns=samp_names)
GT.insert(0, 'AlleleID', LOCI.loc[keep_locus, 'AlleleID'].values)
GT.to_csv(os.path.join(DATA, 'genotypes_filtered_4002x222.csv'), index=False)

log('writing data/loci_filtered_4002.csv')
LOCI.loc[keep_locus].to_csv(os.path.join(DATA, 'loci_filtered_4002.csv'), index=False)

log('writing data/samples_analysed_222.csv')
S = pd.DataFrame({
    'sample': samp_names,
    'dart_order': inc['order'].values,
    'id_in_report': inc['id_report'].values,
    'report_column_index': keep_cols,
    'population': inc['population'].values,
    'region': inc['region'].values,
    'sex': inc['sex'].fillna('').values,
})
S['sex'] = S['sex'].replace({'': 'unknown'})
S.to_csv(os.path.join(DATA, 'samples_analysed_222.csv'), index=False)

log('writing data/ki_individuals_metadata.csv')
KI_META = [
    ('KI1', 'Male',   'Juvenile'),
    ('KI2', 'Female', 'Juvenile'),
    ('KI3', 'Female', 'Adult'),
    ('KI4', 'Male',   'Adult'),
    ('KI5', 'Male',   'Sub-adult'),
    ('KI6', 'Female', 'Adult'),
    ('KI7', 'Male',   'Juvenile'),
    ('KI8', 'Female', 'Juvenile'),
]
KIDF = pd.DataFrame(KI_META, columns=['sample', 'sex', 'age_class'])
KIDF['capture_month'] = '2021-05'
KIDF['capture_river'] = 'Rocky River'
KIDF['capture_area'] = 'Flinders Chase National Park, Kangaroo Island, South Australia'
KIDF['dart_order'] = 'DPla24-9425'
KIDF['coordinate_note'] = ('individual-to-site assignment is not released; the platypus is '
                           'listed Endangered in South Australia')
KIDF.to_csv(os.path.join(DATA, 'ki_individuals_metadata.csv'), index=False)

log('writing data/diversity_by_population.csv, data/diversity_by_individual.csv')
shutil.copy2(os.path.join(OUT, 'diversity_by_population.csv'),
             os.path.join(DATA, 'diversity_by_population.csv'))
shutil.copy2(os.path.join(OUT, 'diversity_by_individual.csv'),
             os.path.join(DATA, 'diversity_by_individual.csv'))

# ---------------------------------------------------------------- results/
log('writing results/')
for fn in ['allelic_richness.json', 'assertions.json', 'filter_cascade.json',
           'fst_pairwise_ci.json', 'joint_summary.json', 'kinship_summary.json',
           'locus_class_recovery.json', 'maf_sensitivity.json', 'old_vs_new_ratios.json',
           'oldnew_locus_bias.json', 'pca_eigen.json', 'polymorphism_matched_n.json',
           'rarefaction.json', 'sfs.json',
           'fst_pairwise.csv', 'kinship_pairs.csv', 'pca_scores.csv',
           'TableS_old_vs_new.csv', 'scaffold_chrom_crosswalk.csv',
           'RESULTS_SUMMARY.md']:
    shutil.copy2(os.path.join(OUT, fn), os.path.join(RES, fn))

# fst_matrix.csv is written by the pipeline with an unnamed index column; give it
# a real header so pandas.read_csv(path) with no arguments is well defined.
FM = pd.read_csv(os.path.join(OUT, 'fst_matrix.csv'), index_col=0)
FM.index.name = 'population'
FM.reset_index().to_csv(os.path.join(RES, 'fst_matrix.csv'), index=False)

# ---------------------------------------------------------------- figures
log('copying figures/')
FIGSRC = os.path.join(REV, 'figures')
nfig = 0
for fn in sorted(os.listdir(FIGSRC)):
    if fn.startswith('Figure_') and (fn.endswith('.png') or fn.endswith('.pdf')):
        shutil.copy2(os.path.join(FIGSRC, fn), os.path.join(DEP, 'figures', fn))
        nfig += 1
log(f'  {nfig} figure files copied')

# ---------------------------------------------------------------- code/
log('copying code/')
for fn in ['ki_coprocessed_joint.py', 'build_revision_figures_v2.py', 'build_fig1_v2.py']:
    shutil.copy2(os.path.join(WORK, fn), os.path.join(DEP, 'code', fn))
shutil.copy2(os.path.abspath(__file__), os.path.join(DEP, 'code', 'make_deposit_data.py'))

# ---------------------------------------------------------------- data dictionary
log('writing data/data_dictionary.csv')
DD = [
 ('genotypes_filtered_4002x222.csv','AlleleID','DArT allele identifier; the locus key shared with loci_filtered_4002.csv','text'),
 ('genotypes_filtered_4002x222.csv','<sample>','Alternate-allele dosage for that individual at that locus; 222 columns, one per analysed individual (KI1-KI8 plus 214 mainland)','0 = reference homozygote, 1 = heterozygote, 2 = alternate homozygote; no missing values'),
 ('loci_filtered_4002.csv','AlleleID','DArT allele identifier','text'),
 ('loci_filtered_4002.csv','CloneID','DArT clone (tag) identifier; one SNP per CloneID is retained','integer'),
 ('loci_filtered_4002.csv','SNP','SNP position and base change within the tag','text, e.g. 20:C>T'),
 ('loci_filtered_4002.csv','SnpPosition','Position of the SNP within the 69 bp tag','base pairs, 0-based'),
 ('loci_filtered_4002.csv','scaffold','PTTO01 scaffold of the aligned tag (assembly mOrnAna1.pri.v1)','text'),
 ('loci_filtered_4002.csv','scaffold_pos','Alignment position of the tag on the scaffold','base pairs'),
 ('loci_filtered_4002.csv','refseq_majority','Majority RefSeq chromosome accession for the scaffold, from the AlleleID-anchored crosswalk','text, e.g. NC_041735.1'),
 ('loci_filtered_4002.csv','chrom_label','Chromosome label from the crosswalk','1-21, X1-X5, Y1-Y5, unplaced, unanchored'),
 ('loci_filtered_4002.csv','class','Chromosome class from the crosswalk','autosome / sex / unplaced / unanchored'),
 ('loci_filtered_4002.csv','RepAvg','DArT average reproducibility of the locus','proportion 0-1'),
 ('loci_filtered_4002.csv','AlnCnt','Number of reference alignments reported by DArT','count'),
 ('loci_filtered_4002.csv','AlnEvalue','BLAST E-value of the tag alignment','dimensionless'),
 ('loci_filtered_4002.csv','CallRate_report','Call rate reported by DArT across all 376 report samples','proportion 0-1'),
 ('loci_filtered_4002.csv','CallRate_222','Call rate recomputed across the 222 analysed individuals','proportion 0-1'),
 ('loci_filtered_4002.csv','pass_rep','Passed the RepAvg >= 1.0 filter','1 = pass, 0 = fail'),
 ('loci_filtered_4002.csv','pass_aln','Passed the AlnCnt >= 1 filter','1 = pass, 0 = fail'),
 ('loci_filtered_4002.csv','pass_blast','Passed the AlnEvalue <= 1e-20 filter','1 = pass, 0 = fail'),
 ('loci_filtered_4002.csv','pass_autosome','Assigned to an autosome by the crosswalk','1 = pass, 0 = fail'),
 ('loci_filtered_4002.csv','pass_callrate','Call rate == 1.0 over the 222 analysed individuals','1 = pass, 0 = fail'),
 ('loci_filtered_4002.csv','pass_hwe','Did not fail the per-population Hardy-Weinberg test in two or more of the six populations with n >= 19','1 = pass, 0 = fail'),
 ('loci_filtered_4002.csv','pass_secondary','Survived the one-SNP-per-CloneID rule; identical to retained','1 = pass, 0 = fail'),
 ('loci_filtered_4002.csv','retained','In the final analysis set; 1 for every row of this file','1'),
 ('loci_filtered_4002.csv','retained_by_old_cloneid_rule','Would also have been kept by the superseded CloneID-matching rule of the submitted version','1 = yes, 0 = no'),
 ('loci_filtered_4002.csv','p_joint','Alternate-allele frequency across all 222 individuals','proportion 0-1'),
 ('loci_filtered_4002.csv','maf_joint','Minor-allele frequency across all 222 individuals','proportion 0-0.5'),
 ('loci_filtered_4002.csv','p_KI','Alternate-allele frequency in the eight KI individuals','proportion 0-1'),
 ('loci_filtered_4002.csv','p_mainland','Alternate-allele frequency in the 214 mainland individuals','proportion 0-1'),
 ('loci_filtered_4002.csv','poly_KI','Polymorphic within KI','1 = yes, 0 = no'),
 ('loci_filtered_4002.csv','poly_mainland','Polymorphic within the mainland panel','1 = yes, 0 = no'),
 ('loci_filtered_4002.csv','sexlinked_data_driven_flag','Flagged by the data-driven male/female heterozygosity test; reported, never filtered on','1 = flagged, 0 = not'),
 ('samples_analysed_222.csv','sample','Individual identifier used in every other table','KI1-KI8, or the Mijangos et al. (2022) mainland identifier'),
 ('samples_analysed_222.csv','dart_order','DArT order the individual was sequenced in','DPla19-4171 (mainland) or DPla24-9425 (KI)'),
 ('samples_analysed_222.csv','id_in_report','Identifier as spelled in the DArT report header','text; differs from sample for V15a/V15b'),
 ('samples_analysed_222.csv','report_column_index','0-based genotype-column index in the co-processed report; samples are selected by index because seven identifiers occur in both orders as different animals','integer 0-375'),
 ('samples_analysed_222.csv','population','Population after the Mijangos et al. (2022) recoding','one of 11 labels'),
 ('samples_analysed_222.csv','region','River region','KI, SNOWY, U_MURRAY, BORDER'),
 ('samples_analysed_222.csv','sex','Field-recorded sex','M, F, unknown'),
 ('ki_individuals_metadata.csv','sample','KI individual identifier','KI1-KI8'),
 ('ki_individuals_metadata.csv','sex','Field-recorded sex at capture','Male / Female'),
 ('ki_individuals_metadata.csv','age_class','Field-assigned age class at capture','Juvenile / Sub-adult / Adult'),
 ('ki_individuals_metadata.csv','capture_month','Month of capture','YYYY-MM'),
 ('ki_individuals_metadata.csv','capture_river','River of capture','text'),
 ('ki_individuals_metadata.csv','capture_area','Protected area and region of capture','text'),
 ('ki_individuals_metadata.csv','dart_order','DArT order','DPla24-9425'),
 ('ki_individuals_metadata.csv','coordinate_note','Why per-individual coordinates are withheld','text'),
 ('diversity_by_population.csv','population','Population label','one of 11 labels'),
 ('diversity_by_population.csv','region','River region','KI, SNOWY, U_MURRAY, BORDER'),
 ('diversity_by_population.csv','n','Individuals in the population','count'),
 ('diversity_by_population.csv','n_loci','Loci used','4002 for every row'),
 ('diversity_by_population.csv','Ho','Observed heterozygosity','proportion 0-1'),
 ('diversity_by_population.csv','Ho_SE','Standard error of the mean of Ho across loci','proportion'),
 ('diversity_by_population.csv','He_unbiased','Unbiased Nei (1987) gene diversity','proportion 0-1'),
 ('diversity_by_population.csv','He_SE','Standard error of the mean of He across loci','proportion'),
 ('diversity_by_population.csv','He_hw','Biased (Hardy-Weinberg) expected heterozygosity, for reference only','proportion 0-1'),
 ('diversity_by_population.csv','Fis_WC','Weir and Cockerham (1984) F_IS as a ratio of sums over loci','dimensionless'),
 ('diversity_by_population.csv','Fis_jack_SE','Delete-one-locus jackknife standard error of Fis_WC','dimensionless'),
 ('diversity_by_population.csv','Fis_boot_lo','Lower bound, 1,000-replicate bootstrap-over-loci 95% CI for F_IS','dimensionless'),
 ('diversity_by_population.csv','Fis_boot_hi','Upper bound of the same interval','dimensionless'),
 ('diversity_by_population.csv','n_polymorphic','Loci polymorphic within the population','count of 4002'),
 ('diversity_by_population.csv','prop_polymorphic','n_polymorphic / 4002','proportion 0-1'),
 ('diversity_by_population.csv','Ar_g16','Rarefied allelic richness at g = 16 gene copies (El Mousadik and Petit 1996)','alleles per locus; blank for the two n = 4 populations'),
 ('diversity_by_population.csv','Ar_g16_SE','Standard error of the mean of Ar across loci','alleles per locus'),
 ('diversity_by_population.csv','mean_folded_MAF','Mean folded minor-allele frequency over loci polymorphic in the population','proportion 0-0.5'),
 ('diversity_by_individual.csv','sample','Individual identifier','text'),
 ('diversity_by_individual.csv','population','Population label','one of 11 labels'),
 ('diversity_by_individual.csv','n_loci_called','Loci scored in the individual','4002 for every row'),
 ('diversity_by_individual.csv','Ho_ind','Proportion of scored loci at which the individual is heterozygous','proportion 0-1'),
 ('diversity_by_individual.csv','F_hom_ownpop','Excess-homozygosity F (Yang et al. 2010 F-hat-2) in the individual own-population allele-frequency frame','dimensionless'),
 ('diversity_by_individual.csv','F_uni_ownpop','Correlation of uniting gametes (Yang et al. 2010 F-hat-3), own-population frame','dimensionless'),
 ('diversity_by_individual.csv','F_grm_ownpop','Variance-standardised GRM diagonal minus 1 (F-hat-1), own-population frame','dimensionless'),
 ('diversity_by_individual.csv','F_uni_ownpop_lo95','Lower bound, 1,000-replicate bootstrap-over-loci 95% CI for F_uni_ownpop','dimensionless'),
 ('diversity_by_individual.csv','F_uni_ownpop_hi95','Upper bound of the same interval','dimensionless'),
 ('diversity_by_individual.csv','n_loci_F_ownpop','Loci polymorphic in the own-population frame and used for its F estimators','count'),
 ('diversity_by_individual.csv','F_hom_mainland','Excess-homozygosity F in the pooled 214-sample mainland frame','dimensionless'),
 ('diversity_by_individual.csv','F_uni_mainland','Correlation of uniting gametes, mainland frame','dimensionless'),
 ('diversity_by_individual.csv','F_grm_mainland','GRM-diagonal F, mainland frame; interpret as a rank ordering only','dimensionless'),
 ('diversity_by_individual.csv','F_uni_mainland_lo95','Lower bound, bootstrap 95% CI for F_uni_mainland','dimensionless'),
 ('diversity_by_individual.csv','F_uni_mainland_hi95','Upper bound of the same interval','dimensionless'),
 ('diversity_by_individual.csv','n_loci_F_mainland','Loci polymorphic in the mainland frame and used for its F estimators','count'),
 ('diversity_by_individual.csv','mean_kinship_within_pop','Mean Yang et al. (2010) relatedness of the individual to the others in its population','dimensionless'),
 ('diversity_by_individual.csv','max_kinship_within_pop','Maximum of the same','dimensionless'),
 ('../results/fst_pairwise.csv','pop_i, pop_j','The two populations compared','labels'),
 ('../results/fst_pairwise.csv','n_i, n_j','Individuals in each population','count'),
 ('../results/fst_pairwise.csv','n_loci','Loci used','4002'),
 ('../results/fst_pairwise.csv','Fst_WC','Weir and Cockerham (1984) F_ST, ratio of sums over loci','dimensionless'),
 ('../results/fst_pairwise.csv','Fst_lo95, Fst_hi95','1,000-replicate bootstrap-over-loci 95% CI','dimensionless'),
 ('../results/fst_pairwise.csv','Gst_Nei','Nei G_ST','dimensionless'),
 ('../results/fst_pairwise.csv','GstPP_Hedrick','Hedrick standardised G"_ST','dimensionless'),
 ('../results/fst_pairwise.csv','Hs, Ht','Within- and total-population gene diversity','proportion 0-1'),
 ('../results/fst_matrix.csv','population','Row population label','labels'),
 ('../results/fst_matrix.csv','<population>','Weir and Cockerham F_ST between the row and column populations; diagonal is 0','dimensionless'),
 ('../results/kinship_pairs.csv','sample_i, sample_j','The two individuals compared','text'),
 ('../results/kinship_pairs.csv','population_i, population_j','Their populations','labels'),
 ('../results/kinship_pairs.csv','within_pop','Both individuals in the same population','1 = yes, 0 = no'),
 ('../results/kinship_pairs.csv','kinship_yang2010','Yang et al. (2010) genomic relationship, computed with within-population allele frequencies','dimensionless'),
 ('../results/kinship_pairs.csv','relatedness_scale','The same value on the relatedness scale (unrelated ~0, full sibs and parent-offspring ~0.5); halve for a kinship coefficient','dimensionless'),
 ('../results/kinship_pairs.csv','n_loci','Loci polymorphic in the frame used for the pair','count'),
 ('../results/pca_scores.csv','sample, population, region','Individual and its grouping','text'),
 ('../results/pca_scores.csv','PC1-PC6','Principal component scores under the Patterson et al. (2006) normalisation of the 4,002 x 222 dosage matrix','dimensionless'),
 ('../results/TableS_old_vs_new.csv','population, n','Population and sample size','text, count'),
 ('../results/TableS_old_vs_new.csv','loci_old, loci_new','Loci in the submitted (1,296) and co-processed (4,002) analyses','count'),
 ('../results/TableS_old_vs_new.csv','Ho_old ... Ar16_pct_change','Submitted versus recomputed Ho, He, Fis, polymorphic loci and Ar(16), with differences and percentage changes; the Fis_estimator_* columns name the estimator each value came from','mixed'),
 ('../results/scaffold_chrom_crosswalk.csv','scaffold','PTTO01 scaffold identifier','text'),
 ('../results/scaffold_chrom_crosswalk.csv','refseq_majority','Majority RefSeq chromosome accession among the anchor loci on the scaffold','text'),
 ('../results/scaffold_chrom_crosswalk.csv','chrom_label','Chromosome label for the scaffold','1-21, X1-X5, Y1-Y5, unplaced, unanchored'),
 ('../results/scaffold_chrom_crosswalk.csv','purity','Proportion of anchor loci on the scaffold agreeing with refseq_majority','proportion 0-1'),
 ('../results/scaffold_chrom_crosswalk.csv','n_anchor_loci','Anchor loci (AlleleIDs shared with the Mijangos-only report) on the scaffold','count'),
 ('../results/scaffold_chrom_crosswalk.csv','class','Chromosome class assigned to the scaffold','autosome / sex / unplaced / unanchored'),
]
pd.DataFrame(DD, columns=['file', 'column', 'definition', 'units_or_codes']).to_csv(
    os.path.join(DATA, 'data_dictionary.csv'), index=False)

log('done')
