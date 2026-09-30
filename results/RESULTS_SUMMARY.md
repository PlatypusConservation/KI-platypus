# KI platypus — co-processed re-analysis: results summary

Generated 2026-09-12T01:17:30 by `ki_coprocessed_joint.py` (seed 42; numpy 2.2.6, pandas 2.3.3).  Runtime 10.6 s.

Source: DArT **co-processed** one-row-per-locus mapping report `Report_DPla24-9425_1_moreOrders_SNP_mapping_2.csv` (22,054 loci x 376 samples, both DArT orders called in a single pipeline run).  **No CloneID merge and no batch-effect call-string filter are applied anywhere in this pipeline.**

> All 21 hard assertions specified in the analysis plan passed; values are in `assertions.json` and in section 13 below.


## 1. Filter cascade (supplementary Table S1)

| # | Filter | Loci passing (standalone) | Loci remaining (cumulative) |
|---|---|---|---|
| 0 | input | 22,054 | **22,054** |
| 1 | RepAvg>=1.0 | 14,795 | **14,795** |
| 2 | AlnCnt>=1 | 20,770 | **13,966** |
| 3 | AlnEvalue<=1e-20 | 20,040 | **13,461** |
| 4 | autosomal (crosswalk) | 19,484 | **11,619** |
| 5 | CallRate==1.0 over 222 | 7,537 | **4,163** |
| 6 | HWE per-population (fail in >=2 of 6 pops, Bonferroni) | 3 removed | **4,160** |
| 7 | one SNP per CloneID | one per CloneID | **4,002** |

**Final locus set: 4,002 SNPs x 222 individuals, zero missing genotypes** (3,473 polymorphic in the joint panel, 529 monomorphic).  This is 3.09x the 1,296 loci of the submitted analysis.

Sex-chromosome exclusion removed **2,570 loci (11.7%)** via a scaffold->chromosome crosswalk built from 12,077 AlleleID anchors shared with the Mijangos-only report (376 of 581 scaffolds anchored, 359 of them >=95% pure for one chromosome; 99.92% locus-level agreement with the old report's direct annotation).

Sensitivities: 
- HWE form — per-population (primary, 3 loci removed) gives 4,002 loci; the joint-sample Bonferroni form quoted in the submitted Methods removes 288 and gives 3,727.
- Unanchored scaffolds — dropping the 347 loci on scaffolds with no crosswalk anchor gives 3,982 loci and shifts no population H_E by more than 0.00031.
- Call rate — 1.00: 4,002 loci (3,473 polymorphic); 0.99: 5,314 loci (4,303 polymorphic); 0.98: 5,855 loci (4,555 polymorphic); 0.95: 6,985 loci (4,989 polymorphic); 0.90: 8,026 loci (5,475 polymorphic).
- Residual sex-linkage — 355 loci trip a data-driven male/female heterozygosity test; 8 of them survive into the final set. Reported, not filtered on.


## 2. Per-population diversity (main-text Table 1)

| Population | Region | n | H_O (SE) | H_E (SE) | F_IS (W&C) | jackknife SE | 95% bootstrap CI | Polymorphic loci | % poly | A_r (g=16) |
|---|---|---|---|---|---|---|---|---|---|---|
| **KI** | KI | 8 | 0.0989 (0.0031) | 0.0954 (0.0027) | -0.0385 | 0.0130 | [-0.0627, -0.0131] | 1,105 | 27.6% | 1.2761 |
| EUCUMBENE_ABOVE | SNOWY | 4 | 0.1140 (0.0033) | 0.1135 (0.0030) | -0.0048 | 0.0142 | [-0.0317, +0.0220] | 1,190 | 29.7% | — |
| EUCUMBENE_BELOW | SNOWY | 20 | 0.1090 (0.0028) | 0.1078 (0.0026) | -0.0114 | 0.0070 | [-0.0252, +0.0025] | 1,648 | 41.2% | 1.3436 |
| MITTA_ABOVE | U_MURRAY | 13 | 0.1071 (0.0027) | 0.1101 (0.0026) | +0.0287 | 0.0083 | [+0.0123, +0.0441] | 1,681 | 42.0% | 1.3687 |
| MITTA_BELOW | U_MURRAY | 4 | 0.1174 (0.0033) | 0.1199 (0.0030) | +0.0239 | 0.0145 | [-0.0051, +0.0526] | 1,254 | 31.3% | — |
| OVENS | U_MURRAY | 19 | 0.1123 (0.0026) | 0.1133 (0.0026) | +0.0093 | 0.0066 | [-0.0033, +0.0227] | 1,859 | 46.5% | 1.3786 |
| SEVERN_ABOVE | BORDER | 23 | 0.1058 (0.0027) | 0.1051 (0.0026) | -0.0069 | 0.0062 | [-0.0191, +0.0049] | 1,572 | 39.3% | 1.3270 |
| SEVERN_BELOW | BORDER | 17 | 0.1022 (0.0026) | 0.1030 (0.0026) | +0.0079 | 0.0074 | [-0.0069, +0.0237] | 1,516 | 37.9% | 1.3268 |
| SNOWY | SNOWY | 56 | 0.1102 (0.0025) | 0.1111 (0.0025) | +0.0083 | 0.0042 | [+0.0001, +0.0162] | 2,183 | 54.5% | 1.3758 |
| TENTERFIELD | BORDER | 39 | 0.1055 (0.0025) | 0.1075 (0.0025) | +0.0183 | 0.0049 | [+0.0092, +0.0274] | 1,890 | 47.2% | 1.3523 |
| THREDBO | SNOWY | 19 | 0.1115 (0.0026) | 0.1121 (0.0025) | +0.0051 | 0.0067 | [-0.0083, +0.0180] | 1,839 | 46.0% | 1.3749 |

All 4,002 loci, no missing data. H_E is the unbiased Nei (1987) estimator; F_IS is Weir & Cockerham (1984) as a ratio of sums over loci with delete-one-locus jackknife SE and a 1,000-replicate bootstrap over loci.  A_r is El Mousadik & Petit (1996) rarefied allelic richness at g = 16 gene copies; the two n = 4 populations cannot be standardised to g = 16.

- KI H_E = **0.0954** = **86.5%** of the mainland mean (0.1103) and **84.2%** of Ovens (0.1133).
- KI H_O = 0.0989 = 90.3% of the mainland mean.
- KI A_r(16) = 1.2761 = 94.1% of the mainland mean over the 8 populations that can be standardised.
- KI is the lowest of the eleven populations for H_E (True), H_O (True) and A_r(16) (True).

> **Absolute H_E is not comparable to Mijangos et al. (2022) Table 1 and the submitted cross-check must be deleted.**  The mainland mean here is 0.1103 against the published 0.140 (-21.2%). The co-processed panel is ascertained across all 376 samples and retains loci that are monomorphic within the Mijangos 218, which neither the Mijangos-only report nor the CloneID intersection ever contained. This is expected and is not a bug; every KI-versus-mainland claim below is made in relative terms.


## 3. Reviewer 1 — what the CloneID-matching approach lost, and in which direction

Old rule reconstructed exactly: keep a CloneID iff it is present in both reports **and** the two SNP call strings match (11,451 shared CloneIDs, 10,336 with matching call strings).  Applied to the 4,002 final loci:

| Locus class | n loci | % of final set | Polymorphic in KI (n=8) | Polymorphic in mainland (n=214) | Mean joint MAF | Mean H_E (222) |
|---|---|---|---|---|---|---|
| Retained by the old CloneID rule | 2,953 | 73.8% | 910 (30.8%) | 2,938 (99.5%) | 0.0960 | 0.1458 |
| **Lost by the old CloneID rule** | **1,049** | **26.2%** | 195 (18.6%) | 476 (45.4%) | 0.0397 | 0.0581 |

The lost loci are overwhelmingly the low-frequency, near-fixed ones: 49.1% of them are monomorphic across all 222 individuals versus 0.5% of the retained loci.

**H_E computed on each subset (direction and size of the inflation):**

| Population | All final loci | Retained-by-old subset | Lost-by-old subset | Inflation from using only the old subset |
|---|---|---|---|---|
| KI | 0.0954 | 0.1068 | 0.0633 | +12.0% |
| EUCUMBENE_ABOVE | 0.1135 | 0.1341 | 0.0558 | +18.1% |
| EUCUMBENE_BELOW | 0.1078 | 0.1285 | 0.0497 | +19.2% |
| MITTA_ABOVE | 0.1101 | 0.1315 | 0.0498 | +19.5% |
| MITTA_BELOW | 0.1199 | 0.1433 | 0.0541 | +19.5% |
| OVENS | 0.1133 | 0.1361 | 0.0491 | +20.1% |
| SEVERN_ABOVE | 0.1051 | 0.1257 | 0.0472 | +19.6% |
| SEVERN_BELOW | 0.1030 | 0.1222 | 0.0487 | +18.7% |
| SNOWY | 0.1111 | 0.1325 | 0.0508 | +19.3% |
| TENTERFIELD | 0.1075 | 0.1286 | 0.0479 | +19.7% |
| THREDBO | 0.1121 | 0.1340 | 0.0505 | +19.5% |

**Submitted versus recomputed, per population (supplementary Table S2; full file `TableS_old_vs_new.csv`):**

| Population | n | H_O submitted | H_O new | H_E submitted | H_E new | ΔH_E % | F_IS submitted | F_IS new (W&C) | Polymorphic submitted | Polymorphic new | A_r(16) submitted | A_r(16) new |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| KI | 8 | 0.1337 | 0.0989 | 0.1309 | 0.0954 | -27.1% | -0.0231 | -0.0385 | 472 | 1,105 | 1.3642 | 1.2761 |
| EUCUMBENE_ABOVE | 4 | 0.1609 | 0.1140 | 0.1548 | 0.1135 | -26.6% | -0.0305 | -0.0048 | — | 1,190 | — | — |
| EUCUMBENE_BELOW | 20 | 0.1498 | 0.1090 | 0.1473 | 0.1078 | -26.8% | -0.0099 | -0.0114 | 707 | 1,648 | 1.4587 | 1.3436 |
| MITTA_ABOVE | 13 | 0.1455 | 0.1071 | 0.1505 | 0.1101 | -26.8% | 0.0238 | +0.0287 | — | 1,681 | 1.4857 | 1.3687 |
| MITTA_BELOW | 4 | 0.1561 | 0.1174 | 0.1627 | 0.1199 | -26.3% | 0.0323 | +0.0239 | — | 1,254 | — | — |
| OVENS | 19 | 0.1563 | 0.1123 | 0.1563 | 0.1133 | -27.5% | -0.0038 | +0.0093 | 788 | 1,859 | 1.5073 | 1.3786 |
| SEVERN_ABOVE | 23 | 0.1414 | 0.1058 | 0.1423 | 0.1051 | -26.1% | -0.0009 | -0.0069 | 683 | 1,572 | 1.4392 | 1.3270 |
| SEVERN_BELOW | 17 | 0.1442 | 0.1022 | 0.1440 | 0.1030 | -28.5% | -0.0028 | +0.0079 | — | 1,516 | 1.4458 | 1.3268 |
| SNOWY | 56 | 0.1499 | 0.1102 | 0.1512 | 0.1111 | -26.5% | 0.0075 | +0.0083 | 885 | 2,183 | 1.4968 | 1.3758 |
| TENTERFIELD | 39 | 0.1457 | 0.1055 | 0.1488 | 0.1075 | -27.8% | 0.0158 | +0.0183 | 785 | 1,890 | 1.4742 | 1.3523 |
| THREDBO | 19 | 0.1493 | 0.1115 | 0.1511 | 0.1121 | -25.8% | 0.0085 | +0.0051 | 767 | 1,839 | 1.4898 | 1.3749 |

- Locus count 1,296 -> 4,002 (3.09x).
- Absolute H_E falls by 26.9% on average for the mainland populations and by 27.1% for KI.
- **The relative conclusion survives:** KI:mainland-mean H_E 86.7% -> 86.5%; KI:Ovens 83.7% -> 84.2%; KI remains the lowest of the eleven populations.


## 4. Locus-class recovery — Reviewer 1's technical objection

| Locus class | Input panel (22,054) | Final filtered set | Of the final set, would have been dropped by the CloneID rule |
|---|---|---|---|
| Polymorphic within DPla19-4171 (Mijangos order) | 16,773 | 3,420 | 480 (14.0%) |
| Polymorphic within DPla24-9425 (2024 order) | 20,199 | 3,704 | 1,030 (27.8%) |
| Monomorphic in DPla19-4171, polymorphic in DPla24-9425 | 5,222 | 581 | 568 (97.8%) |
| Monomorphic in DPla24-9425, polymorphic in DPla19-4171 | 1,796 | 297 | 18 (6.1%) |
| Monomorphic in both orders | 59 | 1 | 1 (100.0%) |
| Monomorphic in the mainland 214, polymorphic in KI | 566 | 57 | 56 (98.2%) |
| Polymorphic in the mainland 214, fixed in KI | 11,453 | 2,366 | 337 (14.2%) |
| Fixed allelic difference, KI vs mainland | 86 | 2 | 2 (100.0%) |

Because both DArT orders were called in one pipeline run, loci that are monomorphic within one order and loci that are fixed differences between KI and the mainland are present in the matrix and are retained. The CloneID rule discarded a large share of exactly these classes.


## 5. Pairwise differentiation — F_ST (NEW)

| Comparison | F_ST (W&C) | 95% bootstrap CI | G_ST (Nei) | G"_ST (Hedrick) |
|---|---|---|---|---|
| KI vs OVENS | 0.1468 | [0.1359, 0.1601] | 0.0826 | 0.1703 |
| KI vs MITTA_BELOW | 0.1520 | [0.1372, 0.1674] | 0.0775 | 0.1611 |
| KI vs MITTA_ABOVE | 0.1558 | [0.1435, 0.1685] | 0.0863 | 0.1772 |
| KI vs THREDBO | 0.1599 | [0.1473, 0.1745] | 0.0902 | 0.1847 |
| KI vs SNOWY | 0.1618 | [0.1489, 0.1760] | 0.0935 | 0.1908 |
| KI vs EUCUMBENE_BELOW | 0.1820 | [0.1684, 0.1957] | 0.1030 | 0.2079 |
| KI vs EUCUMBENE_ABOVE | 0.1884 | [0.1726, 0.2056] | 0.1000 | 0.2030 |
| KI vs TENTERFIELD | 0.2404 | [0.2256, 0.2558] | 0.1419 | 0.2767 |
| KI vs SEVERN_ABOVE | 0.2602 | [0.2440, 0.2758] | 0.1530 | 0.2950 |
| KI vs SEVERN_BELOW | 0.2637 | [0.2485, 0.2797] | 0.1540 | 0.2963 |

- KI is **most similar to OVENS** (F_ST = 0.1468) and most divergent from SEVERN_BELOW (F_ST = 0.2637).
- Mean KI-vs-mainland F_ST = 0.1911; mean mainland-vs-mainland F_ST = 0.1115 (range 0.0207–0.2015).

The full 11 x 11 matrix is in `fst_matrix.csv` and `fst_pairwise_ci.json`.

| | KI | EUCUMBENE_ABOVE | EUCUMBENE_BELOW | MITTA_ABOVE | MITTA_BELOW | OVENS | SEVERN_ABOVE | SEVERN_BELOW | SNOWY | TENTERFIELD | THREDBO |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **KI** | — | 0.188 | 0.182 | 0.156 | 0.152 | 0.147 | 0.260 | 0.264 | 0.162 | 0.240 | 0.160 |
| **EUCUMBENE_ABOVE** | 0.188 | — | 0.053 | 0.070 | 0.052 | 0.072 | 0.191 | 0.197 | 0.039 | 0.178 | 0.037 |
| **EUCUMBENE_BELOW** | 0.182 | 0.053 | — | 0.087 | 0.075 | 0.088 | 0.198 | 0.202 | 0.037 | 0.184 | 0.027 |
| **MITTA_ABOVE** | 0.156 | 0.070 | 0.087 | — | 0.021 | 0.045 | 0.172 | 0.173 | 0.069 | 0.162 | 0.069 |
| **MITTA_BELOW** | 0.152 | 0.052 | 0.075 | 0.021 | — | 0.030 | 0.172 | 0.172 | 0.056 | 0.156 | 0.051 |
| **OVENS** | 0.147 | 0.072 | 0.088 | 0.045 | 0.030 | — | 0.169 | 0.170 | 0.071 | 0.159 | 0.066 |
| **SEVERN_ABOVE** | 0.260 | 0.191 | 0.198 | 0.172 | 0.172 | 0.169 | — | 0.046 | 0.179 | 0.063 | 0.180 |
| **SEVERN_BELOW** | 0.264 | 0.197 | 0.202 | 0.173 | 0.172 | 0.170 | 0.046 | — | 0.179 | 0.060 | 0.183 |
| **SNOWY** | 0.162 | 0.039 | 0.037 | 0.069 | 0.056 | 0.071 | 0.179 | 0.179 | — | 0.169 | 0.021 |
| **TENTERFIELD** | 0.240 | 0.178 | 0.184 | 0.162 | 0.156 | 0.159 | 0.063 | 0.060 | 0.169 | — | 0.169 |
| **THREDBO** | 0.160 | 0.037 | 0.027 | 0.069 | 0.051 | 0.066 | 0.180 | 0.183 | 0.021 | 0.169 | — |

## 6. PCA — honest answer to the reviewer's overlap objection

Patterson/Price normalisation, SVD on 4,002 loci x 222 individuals, zero missing data. PC1 7.50%, PC2 2.86%, PC3 2.00%, PC4 1.77% of variance.

| Axis | % variance | KI range | Mainland range | Fully separated? | KI individuals inside the mainland range |
|---|---|---|---|---|---|
| PC1 | 7.50% | [+9.48, +12.00] | [-25.68, +15.42] | no | 8 of 8 (KI7, KI8, KI1, KI2, KI3, KI4, KI5, KI6) |
| PC2 | 2.86% | [+44.28, +51.76] | [-8.57, +10.54] | **yes** | 0 of 8 (none) |
| PC3 | 2.00% | [+16.75, +21.95] | [-24.65, +8.07] | **yes** | 0 of 8 (none) |
| PC4 | 1.77% | [-4.52, -1.15] | [-19.30, +17.83] | no | 8 of 8 (KI7, KI8, KI1, KI2, KI3, KI4, KI5, KI6) |

**Statement for the manuscript.** KI is completely separated from all 214 mainland individuals on PC2, PC3 (on PC2 the KI minimum is +44.28 against a mainland maximum of +10.54, a gap of 33.74 units), and is **not** separable on PC1, PC4, where every KI individual falls inside the mainland spread. On the PC1 x PC2 plane the minimum KI-to-mainland Euclidean distance is 33.81 against a maximum within-KI distance of 7.54, so the eight KI animals form a discrete cluster in that plane with no mainland individual closer to any KI individual than the KI animals are to each other.

Because PC1 is dominated by the mainland Border-versus-Snowy contrast, the axis that isolates KI is PC2, not PC1. The figure must therefore be drawn and described as PC1 x PC2 with KI separated on PC2; the submitted claim that KI forms "a discrete cluster well removed" from the mainland should be replaced by this explicit, axis-by-axis statement. No KI individual falls inside the mainland spread on PC2, PC3.

Population means on the first four axes:

| Population | PC1 | PC2 | PC3 | PC4 |
|---|---|---|---|---|
| KI | +10.96 | +47.75 | +18.97 | -2.65 |
| EUCUMBENE_ABOVE | +12.90 | -3.46 | +0.05 | +0.02 |
| EUCUMBENE_BELOW | +13.90 | -6.55 | +5.49 | +0.11 |
| MITTA_ABOVE | +9.71 | +6.32 | -16.10 | +1.81 |
| MITTA_BELOW | +10.40 | +6.10 | -15.59 | +2.08 |
| OVENS | +9.51 | +7.98 | -19.05 | +2.02 |
| SEVERN_ABOVE | -22.83 | -0.95 | -1.06 | -14.56 |
| SEVERN_BELOW | -22.67 | -0.86 | -0.69 | -10.83 |
| SNOWY | +13.56 | -6.25 | +4.40 | -0.76 |
| TENTERFIELD | -22.23 | -0.51 | +2.43 | +13.34 |
| THREDBO | +13.26 | -4.67 | +3.51 | -0.50 |

## 7. Individual-level heterozygosity and inbreeding (Reviewer 2)

| Sample | H_O | F_HOM (KI frame) | F_UNI (KI frame) | F_GRM (KI frame) | F_UNI 95% CI (KI frame) | F_HOM (mainland frame) | F_UNI (mainland frame) | F_GRM (mainland frame) |
|---|---|---|---|---|---|---|---|---|
| KI7 | 0.1029 | -0.0787 | -0.1483 | -0.2427 | [-0.1905, -0.1024] | +0.2082 | +0.4011 | +1.4264 |
| KI8 | 0.0865 | +0.0941 | -0.0091 | -0.0667 | [-0.0609, +0.0402] | +0.3449 | +0.8279 | +1.8512 |
| KI1 | 0.1024 | -0.0735 | -0.1413 | -0.1869 | [-0.1858, -0.0978] | +0.2143 | +0.5494 | +1.7271 |
| KI2 | 0.0967 | -0.0132 | -0.0734 | -0.1009 | [-0.1238, -0.0208] | +0.2531 | +0.3344 | +1.0564 |
| KI3 | 0.1004 | -0.0525 | -0.1284 | -0.1755 | [-0.1741, -0.0874] | +0.2347 | +0.2579 | +1.1437 |
| KI4 | 0.0972 | -0.0185 | -0.0257 | +0.2373 | [-0.0843, +0.0339] | +0.2572 | +0.3543 | +1.2784 |
| KI5 | 0.1042 | -0.0918 | -0.1305 | -0.0651 | [-0.1771, -0.0852] | +0.2061 | +0.5612 | +1.5655 |
| KI6 | 0.1004 | -0.0525 | -0.1278 | -0.1840 | [-0.1718, -0.0852] | +0.2245 | +0.3267 | +1.1073 |

Computed on all 4,002 loci (F estimators use the 1,105 loci polymorphic in KI and the 3,414 loci polymorphic in the mainland panel respectively). F_HOM = excess homozygosity; F_UNI = Yang et al. (2010) correlation of uniting gametes (F-hat-3); F_GRM = variance-standardised GRM diagonal minus 1 (F-hat-1). CIs are 1,000-replicate bootstraps over loci.

- Per-individual H_O: KI 0.0865–0.1042 (mean 0.0989); mainland 0.0720–0.1439 (mean 0.1084). Every KI animal sits below the mainland mean and the KI mean is at the 4.7th percentile of the mainland distribution, but the KI and mainland ranges **overlap**: 45 of the 214 mainland animals have a lower H_O than the highest KI animal and 0 KI animal(s) fall below the mainland minimum. The KI deficit is a shift in the whole distribution, not a set of individually exceptional animals — state it that way.
- In the **mainland reference frame** (the frame that puts all 222 animals on a common scale) the KI animals have F_UNI +0.258 to +0.828 (mean +0.452) against a mainland mean of +0.079 and a mainland maximum of +0.370; 4 of 8 exceed the highest mainland value.
- F_GRM (F-hat-1) is the most rare-allele-sensitive of the three estimators and inflates strongly when evaluated in a foreign reference frame; the mainland-frame F_GRM values above should be read as a rank ordering, not as an absolute inbreeding coefficient. F_HOM and F_UNI are the ones to quote.
- In the **KI reference frame** the same animals have F_UNI -0.148 to -0.009. Within-sample frequencies force this mean towards zero by construction, so own-frame values describe variation *among* the eight animals, not their absolute level; the matched n = 8 mainland null for F_UNI is -0.062 (-0.117 to -0.008, 5th–95th percentile over 200 draws from each of the six populations with n >= 19).


## 8. Pairwise kinship among the eight KI individuals (Reviewer 2)

Yang et al. (2010) genomic relationship, on the **relatedness** scale (unrelated ~0, full sibs and parent–offspring ~0.5); halve for the kinship coefficient. KI values use KI allele frequencies.

| | KI7 | KI8 | KI1 | KI2 | KI3 | KI4 | KI5 | KI6 |
|---|---|---|---|---|---|---|---|---|
| **KI7** | — | -0.265 | +0.203 | -0.225 | -0.238 | -0.185 | +0.222 | -0.269 |
| **KI8** | -0.265 | — | -0.271 | +0.089 | -0.095 | -0.218 | -0.294 | +0.121 |
| **KI1** | +0.203 | -0.271 | — | -0.196 | -0.281 | -0.183 | +0.204 | -0.290 |
| **KI2** | -0.225 | +0.089 | -0.196 | — | -0.036 | -0.158 | -0.277 | -0.097 |
| **KI3** | -0.238 | -0.095 | -0.281 | -0.036 | — | -0.098 | -0.286 | +0.209 |
| **KI4** | -0.185 | -0.218 | -0.183 | -0.158 | -0.098 | — | -0.205 | -0.191 |
| **KI5** | +0.222 | -0.294 | +0.204 | -0.277 | -0.286 | -0.205 | — | -0.298 |
| **KI6** | -0.269 | +0.121 | -0.290 | -0.097 | +0.209 | -0.191 | -0.298 | — |

- KI: 28 pairs, mean -0.1288, range -0.2982 to +0.2218.
- Matched n = 8 mainland null (200 draws from each of EUCUMBENE_BELOW, OVENS, SEVERN_ABOVE, SNOWY, TENTERFIELD, THREDBO; 33,600 pairs): mean -0.1340, 95th percentile -0.0546, 99th +0.1540, maximum +0.2982.
- 7 of the 28 KI pairs exceed the matched-n null 95th percentile and 0 exceed its maximum. The empirical p for the single highest KI pair is 0.0031.
- On the full-sample empirical null (all within-population pairs in the six populations with n >= 19, 3,066 pairs): mean -0.0284, 95th +0.0319, 99.9th +0.4286, maximum +0.5301.
- KI pairs above the matched n = 8 mainland 95th percentile (-0.0546): KI7–KI5 (+0.222), KI3–KI6 (+0.209), KI1–KI5 (+0.204), KI7–KI1 (+0.203), KI8–KI6 (+0.121), KI8–KI2 (+0.089), KI2–KI3 (-0.036). Only the four positive values are plausibly informative; the highest three (KI1, KI5 and KI7 in a triangle) are the closest thing to structure in the KI sample.
- **Verdict:** no KI pair reaches a relatedness consistent with first-order kinship (4 of 28 pairs exceed 0.125, the nominal second-order threshold). The claim should be stated against the matched-n null shown above rather than against a theoretical threshold, because within-sample allele frequencies at n = 8 shift the whole distribution downward by construction.


## 9. Proportion of polymorphic loci — replacement metric

The submitted metric ("3,318 of 22,054 input loci, ~15%") is deleted: the denominator was the unfiltered panel, singletons were weighted like intermediate-frequency sites, and n = 8 makes an unobserved alternate allele unremarkable. Three replacement statements, in order:

**(1) On the filtered joint set**, for all eleven populations (denominator 4,002):

| Population | n | Polymorphic loci | % |
|---|---|---|---|
| SNOWY | 56 | 2,183 | 54.5% |
| TENTERFIELD | 39 | 1,890 | 47.2% |
| OVENS | 19 | 1,859 | 46.5% |
| THREDBO | 19 | 1,839 | 46.0% |
| MITTA_ABOVE | 13 | 1,681 | 42.0% |
| EUCUMBENE_BELOW | 20 | 1,648 | 41.2% |
| SEVERN_ABOVE | 23 | 1,572 | 39.3% |
| SEVERN_BELOW | 17 | 1,516 | 37.9% |
| MITTA_BELOW | 4 | 1,254 | 31.3% |
| EUCUMBENE_ABOVE | 4 | 1,190 | 29.7% |
| KI | 8 | 1,105 | 27.6% |

**(2) Sample-size matched at n = 8** (200 draws; this is the number for the abstract):

| Population | Polymorphic loci at n = 8 (mean ± SD) | KI as % of this | KI deficit |
|---|---|---|---|
| EUCUMBENE_BELOW | 1376 ± 31 | 80.3% | 19.7% |
| OVENS | 1515 ± 26 | 73.0% | 27.0% |
| SEVERN_ABOVE | 1304 ± 31 | 84.8% | 15.2% |
| SNOWY | 1501 ± 25 | 73.6% | 26.4% |
| TENTERFIELD | 1409 ± 27 | 78.4% | 21.6% |
| THREDBO | 1496 ± 20 | 73.9% | 26.1% |
| **KI (observed, n = 8)** | **1,105** | — | — |

Against the mean of the six matched mainland samples (1433 loci), KI retains 77.1% — a sample-size-controlled deficit of **22.9%**.

**(3) Rarefied allelic richness at g = 16 gene copies** (the properly standardised metric): KI 1.2761 versus a mainland mean of 1.3560 over the 8 populations that can be standardised = 94.1%. Full curves for g = 2…16 are in `allelic_richness.json`.


## 10. Folded site frequency spectra — KI versus a rarefied mainland

The submitted comparison (KI n = 8 versus Ovens n = 19) is confounded: at n = 8 the smallest non-zero folded minor allele frequency is 1/16 = 0.0625, while at n = 19 it is 1/38 = 0.026, so Ovens can express rare-allele classes KI cannot. The primary comparison here is therefore against Ovens rarefied to n = 8.

| Spectrum | n | Polymorphic loci | Mean folded MAF | Median |
|---|---|---|---|---|
| **A. KI** | 8 | 1,105 | **0.2384** | 0.1875 |
| B. Ovens, full (legacy comparison) | 19 | 1,859 | 0.1635 | 0.1316 |
| **C. OVENS, rarefied to n = 8** | 8 | 1513 ± 28 | 0.1962 ± 0.0025 | — |
| C. SNOWY, rarefied to n = 8 | 8 | 1499 ± 27 | 0.1939 ± 0.0031 | — |
| C. THREDBO, rarefied to n = 8 | 8 | 1500 ± 21 | 0.1960 ± 0.0022 | — |
| C. TENTERFIELD, rarefied to n = 8 | 8 | 1407 ± 25 | 0.2007 ± 0.0031 | — |
| D. KI on the loci the old CloneID rule kept | 8 | 910 | 0.2382 | 0.2500 |

**KI folded spectrum, counts by minor allele count out of 16 gene copies:**

| MAC (of 16) | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| KI count | 221 | 178 | 154 | 128 | 134 | 110 | 106 | 74 |
| KI proportion | 0.200 | 0.161 | 0.139 | 0.116 | 0.121 | 0.100 | 0.096 | 0.067 |
| Ovens n=8 proportion | 0.305 | 0.193 | 0.137 | 0.104 | 0.084 | 0.074 | 0.069 | 0.034 |

- **Primary test (A vs C, matched n = 8):** mean KS D over 200 Ovens subsamples = 0.139 ± 0.010; empirical p (proportion of Ovens n = 8 subsamples with a mean folded MAF at least as high as KI's) = 0/200 (< 0.005). Difference in mean folded MAF = +0.0422.
- **Legacy test (A vs B, unmatched n):** KS D = 0.318, p = 5.38e-62; Mann–Whitney z = 15.6, p = 8.35e-55. Reported only for continuity — **44% of this gap is arithmetic**, not biological (the unmatched gap is +0.0749; the matched gap is +0.0422).
- **Reviewer 1's hypothesis (A vs D):** KS D = 0.009, p = 1; KI mean folded MAF on all final loci 0.2384 versus 0.2382 on the old-rule subset.

**Verdict.** On the co-processed data KI's mean folded MAF is 0.238 over 1105 polymorphic loci (submitted: 0.248 over 472), and the folded spectrum is close to monotonically decreasing rather than intermediate-peaked. Restricting KI to the loci the old CloneID rule would have kept changes its mean folded MAF by only -0.0002 (to 0.238), so the merge itself did NOT create the shape: Reviewer 1's specific mechanism is not supported. What DOES explain roughly half of the published KI-Ovens difference is unmatched sample size: rarefying Ovens to n = 8 raises its mean folded MAF from 0.163 to 0.196, removing 44% of the published gap. A residual excess of intermediate-frequency variants in KI remains (KI 0.238 vs the mainland n = 8 envelope 0.194-0.201; empirical p < 0.005), which is consistent with drift at small effective population size but is not diagnostic of an acute recent bottleneck.


## 11. MAF / MAC sensitivity

Thresholds applied symmetrically across the joint 222-sample panel, never within KI alone.

| Scenario | Loci | KI H_E | Mainland mean H_E | Ovens H_E | KI as % of mainland | KI as % of Ovens | KI lowest? |
|---|---|---|---|---|---|---|---|
| none (primary) | 4,002 | 0.0954 | 0.1103 | 0.1133 | 86.5% | 84.2% | yes |
| MAC>=3 | 3,072 | 0.1225 | 0.1427 | 0.1461 | 85.8% | 83.8% | yes |
| MAF>=0.01 | 2,805 | 0.1319 | 0.1547 | 0.1576 | 85.2% | 83.7% | yes |
| MAF>=0.02 | 2,356 | 0.1466 | 0.1787 | 0.1806 | 82.1% | 81.2% | yes |
| MAF>=0.05 | 1,572 | 0.1899 | 0.2391 | 0.2402 | 79.4% | 79.0% | yes |

Absolute H_E rises steeply as rare variants are removed, but the KI:mainland ratio is stable to within a few percentage points and KI remains the lowest population under every threshold.


## 12. Rarefaction

200 subsamples without replacement at n = 5, 8, 10, 15 and 20 from each of the six populations with n >= 19 (seed 42).

| Population | n (full) | H_E full | H_E at n=8 | % of full | H_O at n=8 | F_IS at n=8 | Polymorphic at n=8 | A_r(16) at n=8 |
|---|---|---|---|---|---|---|---|---|
| EUCUMBENE_BELOW | 20 | 0.1078 | 0.1080 ± 0.0013 | 100.1% | 0.1090 | -0.0107 | 1376 ± 31 | 1.3437 |
| OVENS | 19 | 0.1133 | 0.1134 ± 0.0014 | 100.0% | 0.1122 | +0.0106 | 1515 ± 26 | 1.3785 |
| SEVERN_ABOVE | 23 | 0.1051 | 0.1050 ± 0.0018 | 99.9% | 0.1057 | -0.0068 | 1304 ± 31 | 1.3258 |
| SNOWY | 56 | 0.1111 | 0.1110 ± 0.0013 | 99.9% | 0.1102 | +0.0081 | 1501 ± 25 | 1.3751 |
| TENTERFIELD | 39 | 0.1075 | 0.1074 ± 0.0013 | 100.0% | 0.1055 | +0.0195 | 1409 ± 27 | 1.3521 |
| THREDBO | 19 | 0.1121 | 0.1120 ± 0.0009 | 99.9% | 0.1116 | +0.0036 | 1496 ± 20 | 1.3738 |
| **KI (observed)** | **8** | — | **0.0954** | — | **0.0989** | **-0.0385** | **1,105** | **1.2761** |

H_E is essentially unbiased by sample size (every population is within ~1% of its full-sample value at n = 8), so the KI H_E deficit is not a sample-size artefact. Polymorphic-locus counts and allelic richness are strongly sample-size dependent and must be compared at matched n, as in section 9. Full output including n = 5, 10, 15 and 20 is in `rarefaction.json`.


## 13. Assertions (all must pass for the run to be valid)

| Assertion | Value | Expected | Result |
|---|---|---|---|
| `n_columns_402` | 402 | 402 | PASS |
| `n_data_rows_22054` | 22054 | 22054 | PASS |
| `all_rows_402_fields` | 0 | 0 | PASS |
| `meta_col_names` | True | True | PASS |
| `n_samples_376` | 376 | 376 | PASS |
| `order_counts` | {'DPla19-4171': 218, 'DPla24-9425': 158} | {'DPla19-4171': 218, 'DPla24-9425': 158} | PASS |
| `encoding_orientation_lt_1e-5` | 5e-07 | < 1e-5 | PASS |
| `n_samples_retained_222` | 222 | 222 | PASS |
| `n_KI_8` | 8 | 8 | PASS |
| `population_size_table` | {'EUCUMBENE_BELOW': 20, 'SNOWY': 56, 'THREDBO': 19, 'EUCUMBENE_ABOVE': 4, 'TENTERFIELD'... | {'SNOWY': 56, 'TENTERFIELD': 39, 'SEVERN_ABOVE': 23, 'EUCUMBENE_BEL... | PASS |
| `duplicate_ids_across_orders_detected` | ['T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9'] | ['T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9'] | PASS |
| `no_duplicate_id_kept_twice` | ['T4', 'T5', 'T6', 'T7', 'T8', 'T9'] | each duplicated ID retained at most once (only the DPla19-4171 copy) | PASS |
| `sexset_size_10` | ['NC_041749.1', 'NC_041750.1', 'NC_041751.1', 'NC_041752.1', 'NC_041753.1', 'NW_0216380... | 10 RefSeq accessions (X1-X5, Y1-Y5) | PASS |
| `old_report_annotated_alleles` | 14621 | 14621 | PASS |
| `n_sex_loci_2300_2800` | 2570 | 2300..2800 | PASS |
| `crosswalk_agreement_ge_0.99` | 0.99917 | >= 0.99 | PASS |
| `final_locus_count_2500_6000` | 4002 | 2500..6000 | PASS |
| `no_missing_in_final_matrix` | 0 | 0 | PASS |
| `pca_no_missing` | 0 | 0 | PASS |
| `shared_cloneids_11451` | 11451 | 11451 | PASS |
| `cloneid_callstring_match_10336` | 10336 | 10336 | PASS |

## 14. Output files

- `out/RESULTS_SUMMARY.md` — 29,501 bytes
- `out/TableS_old_vs_new.csv` — 4,085 bytes
- `out/allelic_richness.json` — 4,109 bytes
- `out/assertions.json` — 2,822 bytes
- `out/diversity_by_individual.csv` — 66,377 bytes
- `out/diversity_by_population.csv` — 3,341 bytes
- `out/filter_cascade.json` — 3,211 bytes
- `out/fst_matrix.csv` — 2,467 bytes
- `out/fst_pairwise.csv` — 9,484 bytes
- `out/fst_pairwise_ci.json` — 23,842 bytes
- `out/joint_summary.json` — 93,252 bytes
- `out/kinship_pairs.csv` — 251,275 bytes
- `out/kinship_summary.json` — 7,903 bytes
- `out/loci_final.csv` — 4,962,946 bytes
- `out/locus_class_recovery.json` — 2,524 bytes
- `out/maf_sensitivity.json` — 8,633 bytes
- `out/old_vs_new_ratios.json` — 1,170 bytes
- `out/oldnew_locus_bias.json` — 9,010 bytes
- `out/pca_eigen.json` — 4,423 bytes
- `out/pca_scores.csv` — 30,170 bytes
- `out/polymorphism_matched_n.json` — 2,526 bytes
- `out/rarefaction.json` — 15,039 bytes
- `out/samples_used.csv` — 18,132 bytes
- `out/scaffold_chrom_crosswalk.csv` — 31,697 bytes
- `out/scaffold_chrom_crosswalk_DRAFT.csv` — 16,985 bytes
- `out/sfs.json` — 8,465 bytes
