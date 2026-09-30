# Drifting alone: genome-wide diversity in the isolated, introduced Kangaroo Island platypus

Code and data for:

> Bino, G., Hawke, T., Baring, R. & Gongora, J. **Drifting alone: genome-wide diversity in the
> isolated, introduced Kangaroo Island platypus (*Ornithorhynchus anatinus*)**. *Conservation
> Genetics* (in revision).

Repository: <https://github.com/PlatypusConservation/KI-platypus>
Archive and DOI: Zenodo, `10.5281/zenodo.XXXXXXX` (see **How to cite**).
Corresponding author: Gilad Bino, <gilad.bino@unsw.edu.au>, Centre for Ecosystem Science, UNSW Sydney.

---

## 1. What the study is

The platypus (*Ornithorhynchus anatinus*) is endemic to eastern Australia. The Kangaroo Island (KI)
population in South Australia is the species' only population outside its native range. It descends
from approximately 19 individuals released into the Rocky and Breakneck Rivers between 1928 and 1946
— three from Wynyard, Tasmania and 16 from Healesville, Victoria — and has been closed ever since.
It has since experienced streamflow decline, the 2019–20 megafires that burnt 96% of Flinders Chase
National Park, and a record flood; capture rates have roughly halved.

We genotyped eight KI platypuses captured in the Rocky River in May 2021 with DArTseq (order
DPla24-9425), and analysed them against the 214 published mainland Murray–Darling samples of
Mijangos et al. (2022) (order DPla19-4171). The decisive methodological point is that both DArT
orders were **called together in a single pipeline run**, so loci monomorphic within one order, and
fixed allelic differences between the two sets of samples, are present in the matrix rather than
structurally absent from it. Filtering was then applied once, to the joint matrix, so no locus
removal can treat KI and mainland samples differently.

The working dataset is **4,002 biallelic SNPs scored in all 222 individuals with no missing
genotypes**, with every locus assigned to a sex chromosome excluded (3,945 sit on autosome-anchored
scaffolds, 37 on unplaced scaffolds and 20 on scaffolds the crosswalk could not anchor; none carries
an X or Y label). On it, KI has the lowest expected heterozygosity of the eleven populations
(H_E = 0.0954, 86.5% of the mainland mean), the lowest rarefied allelic richness of the nine
populations standardisable to 16 gene copies, 22.9% fewer polymorphic loci at a matched sample size
of eight, and pairwise F_ST of 0.147–0.264 (mean 0.191) against every mainland population, least
against the Victorian Ovens.

---

## 2. What is in this deposit

```
.
├── README.md                     this file
├── LICENSE                       MIT, for the code
├── CITATION.cff                  machine-readable citation metadata
├── .zenodo.json                  Zenodo deposit metadata
├── code/                         the analysis pipeline, the figure scripts, requirements
├── data/                         the filtered genotype matrix and the analysis tables
├── results/                      every statistic the manuscript quotes
└── figures/                      the seven final figures, PNG and PDF
```

Everything the manuscript reports can be recomputed from `data/` alone. `code/reproduce_from_deposit.py`
does exactly that and checks the answers against `results/` (section 4).

### 2.1 `data/` — file-by-file manifest

| File | Rows | Columns | What it is |
|---|---:|---:|---|
| `genotypes_filtered_4002x222.csv` | 4,002 | 223 | **The analysis dataset.** Alternate-allele dosages (0 = reference homozygote, 1 = heterozygote, 2 = alternate homozygote) for the 4,002 retained loci in all 222 analysed individuals. Column 1 is `AlleleID`; columns 2–223 are one individual each, named `KI1`–`KI8` plus the 214 mainland identifiers. No missing values anywhere. |
| `loci_filtered_4002.csv` | 4,002 | 30 | Per-locus metadata for the same 4,002 loci: CloneID, SNP, scaffold and position, crosswalk chromosome assignment, DArT reproducibility and alignment statistics, each filter's pass flag, allele frequencies in the joint panel, in KI and on the mainland, and whether the locus would also have been kept by the superseded CloneID-matching rule. The `class` column is `autosome` for 3,945 loci, `unplaced` for 37 and `unanchored` for 20; no retained locus is sex-linked, and `sexlinked_data_driven_flag` marks the eight that trip the data-driven male/female heterozygosity test and are reported rather than filtered on. |
| `samples_analysed_222.csv` | 222 | 7 | The 222 analysed individuals: identifier, DArT order, identifier as spelled in the report, genotype-column index in the report, population, river region, field-recorded sex. |
| `ki_individuals_metadata.csv` | 8 | 8 | Field metadata for KI1–KI8: sex and age class at capture, capture month, river and protected area. |
| `diversity_by_population.csv` | 11 | 18 | Manuscript **Table 1**. Per-population H_O, unbiased H_E, F_IS with jackknife SE and bootstrap CI, polymorphic-locus count and proportion, rarefied allelic richness at g = 16, mean folded MAF. |
| `diversity_by_individual.csv` | 222 | 18 | Manuscript **Table 2** and Table S8. Per-individual H_O and the three Yang et al. (2010) inbreeding estimators in two reference frames (own population, pooled mainland) with bootstrap CIs, plus mean and maximum within-population relatedness. |
| `data_dictionary.csv` | 110 | 4 | Definition and units for every column of every data and result table. |

### 2.2 `results/` — every number the manuscript quotes

| File | Rows × cols | What it is |
|---|---|---|
| `RESULTS_SUMMARY.md` | — | **The authoritative numbers file.** Fourteen sections covering the filter cascade, per-population diversity, the CloneID-rule diagnostic, locus-class recovery, F_ST, PCA, individual inbreeding, kinship, the polymorphism metrics, the folded spectra, the MAF sensitivity, the rarefaction, and the 21 hard assertions. Every value in the manuscript traces to a section here. |
| `fst_pairwise.csv` | 55 × 12 | Manuscript **Table 3** and Table S6. All 55 pairwise comparisons: Weir & Cockerham F_ST with bootstrap 95% CI, Nei G_ST, Hedrick G″_ST, H_S and H_T. |
| `fst_matrix.csv` | 11 × 12 | The same F_ST values as a symmetric 11 × 11 matrix; column 1 is `population`. |
| `fst_pairwise_ci.json` | — | The F_ST matrix and per-pair bootstrap output in one machine-readable object. |
| `pca_scores.csv` | 222 × 9 | Figure 2. PC1–PC6 for each individual under the Patterson et al. (2006) normalisation. |
| `pca_eigen.json` | — | Variance explained per axis, singular values, and the axis-by-axis KI-versus-mainland overlap test. |
| `kinship_pairs.csv` | 3,320 × 8 | Table S7 and Figure S2. Yang et al. (2010) relatedness for every within- and between-population pair, including all 28 KI pairs. |
| `kinship_summary.json` | — | The KI kinship distribution against the matched n = 8 and full-sample mainland nulls. |
| `TableS_old_vs_new.csv` | 11 × 21 | Table S3. Submitted (1,296-locus) versus co-processed (4,002-locus) H_O, H_E, F_IS, polymorphic loci and A_r(16), per population, with differences and percentage changes. |
| `old_vs_new_ratios.json` | — | The KI-to-mainland ratios under both analyses, which is what shows the relative conclusion surviving. |
| `oldnew_locus_bias.json` | — | What the superseded CloneID rule discarded, and in which direction it biased H_E. |
| `locus_class_recovery.json` | — | Table S2. Locus classes recovered by co-processing, and the share of each the CloneID rule would have dropped. |
| `filter_cascade.json` | — | Table S1. The filter cascade, plus the call-rate, HWE-form and unanchored-scaffold sensitivities. |
| `scaffold_chrom_crosswalk.csv` | 581 × 6 | Table S9. The scaffold-to-chromosome crosswalk used for sex-chromosome exclusion, with per-scaffold purity and anchor counts. |
| `allelic_richness.json` | — | Rarefied allelic richness curves for g = 2…16 gene copies, per population. |
| `polymorphism_matched_n.json` | — | Polymorphic-locus counts at a matched n = 8, per mainland population, and the KI deficit. |
| `rarefaction.json` | — | Table 4 and Figure S1. H_E, H_O, F_IS, polymorphic loci and A_r at n = 5, 8, 10, 15 and 20 for the six populations with n ≥ 19. |
| `sfs.json` | — | Figures 4 and S3. Folded site-frequency spectra for KI, Ovens at full n and rarefied to n = 8, three further rarefied mainland populations, and the two subsets defined by the old CloneID rule, with the matched and legacy tests. |
| `maf_sensitivity.json` | — | Table S4. Diversity and the KI-to-mainland ratio under MAC ≥ 3 and MAF ≥ 0.01, 0.02 and 0.05. |
| `joint_summary.json` | — | One object holding the whole run: every table above plus the run metadata. |
| `assertions.json` | — | The 21 hard assertions the pipeline requires to pass before it writes anything, with value, expectation and result. |

### 2.3 `figures/` — the seven final figures

Each is present as a 300 dpi PNG and as a vector PDF.

| File stem | Figure |
|---|---|
| `Figure_1_study_area` | Two-panel study area: south-eastern Australia with the KI recipient population, the two founder sources and the mainland reference centroids; and the Rocky River with the four capture sites. |
| `Figure_2_PCA` | PCA, two panels (PC1 × PC2, PC2 × PC3) with marginal distributions and 95% ellipses. |
| `Figure_3_diversity_individual` | Per-population diversity with per-individual values overlaid. |
| `Figure_4_SFS_matched_n` | Folded site-frequency spectrum, KI against the mainland envelope at a matched n = 8. |
| `Figure_S1_rarefaction` | Rarefaction of the six largest mainland populations, three panels. |
| `Figure_S2_kinship` | KI pairwise relatedness against the matched n = 8 mainland null. |
| `Figure_S3_SFS_legacy_vs_matched` | The legacy unmatched-n spectrum comparison against the matched-n one. |

### 2.4 `code/`

| File | What it does |
|---|---|
| `ki_coprocessed_joint.py` | The single analysis pipeline. Runs end-to-end from the raw DArT co-processed report with no manual steps and writes everything in `results/` plus the two diversity tables in `data/`. 21 hard assertions must pass or it raises. Seed 42 throughout. Runtime about 11 s. |
| `build_revision_figures_v2.py` | Regenerates Figures 2, 3, 4, S1, S2 and S3 and the tables file from the pipeline's output. |
| `build_fig1_v2.py` | Regenerates Figure 1. |
| `make_deposit_data.py` | Builds `data/` and `results/` for this deposit: reconstructs the 4,002 × 222 genotype matrix from the raw report using the pipeline's own locus and sample masks, cross-checks the reconstruction against the pipeline's published per-individual H_O and per-population polymorphic-locus counts, and copies the result tables and figures. |
| `reproduce_from_deposit.py` | Recomputes the manuscript's headline statistics from `data/genotypes_filtered_4002x222.csv` alone and checks them against `results/`. Needs nothing outside this deposit. |
| `requirements.txt` | The exact package versions used. |

---

## 3. Data provenance and permissions — please read before reusing

**What is here is ours to release.** The eight KI genotypes (DArT order DPla24-9425) are ours, and
every derived table in `data/` and `results/` is a product of this analysis.

**What is deliberately not here.**

1. **The full co-processed DArT report** — `Report_DPla24-9425_1_moreOrders_SNP_mapping_2.csv`,
   22,054 loci × 376 samples — is **not** included. 218 of its 376 samples belong to Mijangos et al.
   (2022), and redistributing them requires the data owners' permission, which we do not hold. The
   full co-processed report is **available from the corresponding author subject to the data owners'
   permission.** The mainland genotypes themselves are already public in their own archive
   (Zenodo, `10.5281/zenodo.7039778`); please cite Mijangos et al. (2022) for them.
2. **Raw DArTseq FASTQ reads** are not deposited, for the same reason and because the genotype
   matrix here reproduces every analysis.
3. **Per-individual capture coordinates.** The platypus is listed Endangered in South Australia. All
   eight animals were captured in the Rocky River in Flinders Chase National Park in May 2021; the
   individual-to-site assignment is not released. Site-level coordinates rounded to approximately
   1 km are available from the corresponding author, and full-precision coordinates to editors and
   reviewers on request.
4. **The GIS layers behind Figure 1** (Rocky River, survey sites, national coastline) are
   third-party or site-sensitive and are not redistributed here. Figure 1 is supplied as PNG and PDF.

Nothing in this deposit contains the full 22,054-locus report, any FASTQ, or any mainland genotype
beyond the 4,002 retained loci in the 214 individuals that this study analysed.

Animal ethics and permits: Flinders University Animal Welfare Committee approval AEC BIOL4071-3;
Ministerial Exemption ME9903147 (Primary Industries and Resources South Australia); permit U27033-1
to undertake scientific research in Flinders Chase National Park (South Australian Department for
Environment and Water). All individuals were released at their site of capture.

---

## 4. How to reproduce every number and every figure

### 4.1 Environment

The published results were produced with:

| | Version |
|---|---|
| Python | 3.10.12 (GCC 11.4.0, Linux x86-64) |
| numpy | 2.2.6 |
| pandas | 2.3.3 |
| matplotlib | 3.10.9 |
| pyshp | required by `build_fig1_v2.py` only |

No scipy, no scikit-learn, no R, no dartR. The chi-square tail, the two-sample
Kolmogorov–Smirnov test, Mann–Whitney U and the binomial coefficients are implemented from
`math` and `numpy` inside the pipeline, so there is no statistical dependency to version-match.
Every stochastic step (bootstraps, rarefaction draws) uses seed 42, so reruns are bit-identical.

```bash
python3 -m venv venv
source venv/bin/activate            # Windows PowerShell: .\venv\Scripts\Activate.ps1
pip install -r code/requirements.txt
```

### 4.2 Path A — reproduce every number from this deposit alone (recommended, ~20 s)

This is the path a reviewer should take. It needs no DArT report, no network and no credentials.

```bash
python code/reproduce_from_deposit.py
```

The script reads `data/genotypes_filtered_4002x222.csv`, recomputes the estimators from the
published formulae, and checks them against the deposited tables. It prints the largest discrepancy
for each check and exits 0 only if every check passes. It verifies:

* the matrix is 4,002 loci × 222 individuals, every cell 0/1/2, no missing values, `KI1`–`KI8` present;
* per-population H_O, unbiased H_E, Weir & Cockerham F_IS, polymorphic-locus counts, rarefied
  allelic richness A_r(16) and mean folded MAF, for all eleven populations, against
  `data/diversity_by_population.csv` (manuscript Table 1);
* per-individual H_O for all 222 individuals against `data/diversity_by_individual.csv` (Table 2);
* Weir & Cockerham F_ST for all 55 pairs against `results/fst_pairwise.csv` (Table 3, Table S6), and
  that `results/fst_matrix.csv` agrees with it;
* PCA variance explained for PC1–PC10 and the PC1–PC6 scores, against `results/pca_eigen.json` and
  `results/pca_scores.csv`, and that KI is fully separated from all 214 mainland individuals on PC2;
* the KI folded minor-allele-count spectrum and mean folded MAF, and the Ovens full-sample mean
  folded MAF, against `results/sfs.json` (Figures 4, S3);
* the abstract's headline values: KI H_E = 0.0954, 86.5% of the mainland mean, lowest of the eleven
  populations; lowest A_r(16) of the nine standardisable populations; F_ST 0.147–0.264, mean 0.191;
  F_IS = −0.0385.

Expected last line:

```
all checks passed: the deposited result tables are reproduced from the deposited genotype matrix
```

The statistics this path does **not** re-derive are the ones whose inputs are not in the deposit:
the filter cascade itself (it starts from the 22,054-locus report) and the bootstrap and rarefaction
intervals (seeded resampling of the deposited matrix would reproduce them, but the script checks
point estimates only, to keep its runtime short). Both are in `results/` with the seed recorded.

### 4.3 Path B — reproduce the whole pipeline from the raw report

This reproduces `results/` and the figure files bit-for-bit, and additionally reproduces the filter
cascade. It needs files that are not in this deposit:

1. The co-processed report `Report_DPla24-9425_1_moreOrders_SNP_mapping_2.csv` — request it from the
   corresponding author (section 3).
2. The Mijangos et al. (2022) metadata files `ID_Pop_Platypus.csv`, `new_pop_assignments.csv`,
   `new_pop_assignments_regions.csv` and `chrom_platypus.csv`, and the Mijangos-only report
   `Report_DPla19-4171_SNP_mapping_2.csv`, which the pipeline uses in exactly two metadata-only
   places (the scaffold-to-chromosome crosswalk and the old-versus-new locus bias diagnostic).

Then:

```bash
# 1. edit the path block at the top of code/ki_coprocessed_joint.py
#    (HOME, GEN, MIJ, WORK, OUT, FIGD, NEWREP, OLDREP) to point at your copies
python code/ki_coprocessed_joint.py          # ~11 s; writes out/ and out/RESULTS_SUMMARY.md
python code/build_revision_figures_v2.py     # Figures 2, 3, 4, S1, S2, S3 and TABLES.md
python code/build_fig1_v2.py                 # Figure 1; needs pyshp and the GIS layers
```

The paths in all three scripts are absolute and were written for the author's machine. They are the
only thing you need to change; nothing else in the scripts is environment-dependent. The pipeline
asserts 21 conditions — report shape, order composition, genotype-encoding orientation, sample
counts per population, crosswalk agreement, final locus count, zero missing data — and raises rather
than producing output if any fails, so a path or file mismatch fails loudly.

`build_fig1_v2.py` additionally needs pyshp and three shapefiles that are not redistributable
(section 3.4), and writes to a hard-coded output directory. Figure 1 cannot be regenerated from this
deposit; it is supplied as PNG and PDF.

### 4.4 Which file each manuscript item comes from

| Manuscript item | File |
|---|---|
| Abstract, all values | `results/RESULTS_SUMMARY.md` §1, §2, §5, §9, §10 |
| Table 1 (population diversity) | `data/diversity_by_population.csv` |
| Table 2 (KI individuals) | `data/diversity_by_individual.csv` + `data/ki_individuals_metadata.csv` |
| Table 3 (KI vs mainland F_ST) | `results/fst_pairwise.csv` |
| Table 4 (rarefaction to n = 8) | `results/rarefaction.json` |
| Table S1 (filter cascade) | `results/filter_cascade.json` |
| Table S2 (locus classes) | `results/locus_class_recovery.json` |
| Table S3 (submitted vs co-processed) | `results/TableS_old_vs_new.csv` |
| Table S4 (MAF sensitivity) | `results/maf_sensitivity.json` |
| Table S6 (full 11 × 11 F_ST) | `results/fst_matrix.csv`, `results/fst_pairwise.csv` |
| Table S7 (KI kinship) | `results/kinship_pairs.csv`, `results/kinship_summary.json` |
| Table S8 (mainland individuals) | `data/diversity_by_individual.csv` |
| Table S9 (crosswalk validation) | `results/scaffold_chrom_crosswalk.csv` |
| Figures 1–4, S1–S3 | `figures/` |

---

## 5. How to cite

**Cite the paper** for the findings:

> Bino, G., Hawke, T., Baring, R. & Gongora, J. Drifting alone: genome-wide diversity in the
> isolated, introduced Kangaroo Island platypus (*Ornithorhynchus anatinus*). *Conservation
> Genetics* (in revision).

**Cite this deposit** for the code and data:

> Bino, G., Hawke, T., Baring, R. & Gongora, J. (2026). *Code and data for "Drifting alone:
> genome-wide diversity in the isolated, introduced Kangaroo Island platypus (Ornithorhynchus
> anatinus)"* (v1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

Replace `10.5281/zenodo.XXXXXXX` with the DOI Zenodo mints. Use the **concept DOI**, which always
resolves to the newest version — see `PUBLISH.md` §2.4. `CITATION.cff` and `.zenodo.json` carry the
same placeholder and must be updated in the same commit.

**Also cite the mainland dataset**, which is not ours:

> Mijangos, J.L., Gruber, B., Berry, O., Pacioni, C. & Georges, A. (2022). Data and code for
> "Contrasting genetic diversity and structure in the platypus across its range". Zenodo.
> https://doi.org/10.5281/zenodo.7039778

---

## 6. Licence

**Code** — everything in `code/` — is released under the **MIT Licence**; see `LICENSE`.

**Data, result tables and figures** — everything in `data/`, `results/` and `figures/` — are released
under **CC0 1.0 Universal** (public domain dedication), which is the Zenodo and Dryad norm for
research data. You may copy, modify and redistribute them for any purpose without permission.
Attribution is not legally required under CC0; we ask for it as a scholarly courtesy, and section 5
says how.

These terms cover only what is in this deposit. They do not extend to the Mijangos et al. (2022)
mainland genotypes, which carry their own licence in their own archive, or to the GIS layers behind
Figure 1.

---

## 7. References for the methods

El Mousadik, A. & Petit, R.J. (1996) *Theoretical and Applied Genetics* 92, 832–839.
Hurlbert, S.H. (1971) *Ecology* 52, 577–586.
Kilian, A. et al. (2012) *Methods in Molecular Biology* 888, 67–89.
Mijangos, J.L., Gruber, B., Berry, O., Pacioni, C. & Georges, A. (2022) *Molecular Ecology* 31, 5406–5425.
Nei, M. (1987) *Molecular Evolutionary Genetics*. Columbia University Press.
Patterson, N., Price, A.L. & Reich, D. (2006) *PLoS Genetics* 2, e190.
Sansaloni, C. et al. (2011) *BMC Proceedings* 5, P54.
Weir, B.S. & Cockerham, C.C. (1984) *Evolution* 38, 1358–1370.
Yang, J. et al. (2010) *Nature Genetics* 42, 565–569.
Zhou, Y. et al. (2021) *Nature* 592, 756–762.

---

## 8. Publishing this repository

`PUBLISH.md` in this folder is the step-by-step for creating the GitHub repository, pushing, cutting
the `v1.0.0` release, and getting Zenodo to archive it and mint the DOI — including how to reserve
the DOI before publishing so the real DOI can go into the manuscript.
