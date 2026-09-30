#!/usr/bin/env python3
"""
build_revision_figures_v2.py
============================
Regenerates every figure and every table of the revised Kangaroo Island platypus
conservation-genomics manuscript from the co-processed re-analysis outputs in
    <work>/out/
and writes
    <Revision>/figures/Figure_*.png  (300 dpi)  and  Figure_*.pdf
    <work>/TABLES.md

Numbering follows REVISION_PLAN.md Part C2 (binding):
    Figure 1  study area (copied unchanged from Fig1_KI_map_REVISED.*; it CANNOT
              be regenerated here and still carries a known year-label defect —
              see the figure1() docstring)
    Figure 2  PCA, two panels (PC1 x PC2, PC2 x PC3) with marginal strips and
              95% ellipses for KI and for each mainland regional group
    Figure 3  per-population diversity with per-individual values overlaid (NEW)
    Figure 4  folded SFS at matched n = 8
    Figure S1 rarefaction (three panels)
    Figure S2 KI pairwise kinship vs matched n = 8 mainland null (NEW)
    Figure S3 legacy unmatched vs matched-n SFS (NEW)

Everything is derived from out/ except the full-sample (n = 19) Ovens folded
spectrum, which out/ stores only as summary statistics.  That one spectrum is
recomputed here directly from the DArT co-processed report restricted to the
4,002 retained loci listed in out/loci_final.csv, and the recomputation is
asserted against the pipeline's own stored values (n_polymorphic = 1859,
mean folded MAF = 0.1634863) before it is used.

British/Australian spelling; colourblind-safe (Okabe-Ito) palette throughout.
No scipy / sklearn dependency.  Python 3 + numpy + pandas + matplotlib only.
"""

import os, json, shutil, sys
from pathlib import Path

WORK = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(WORK / ".mplcache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Ellipse

OUT      = WORK / "out"
DERIVED  = WORK / "derived"
REVISION = WORK.parent
FIGS     = REVISION / "figures"
GENETICS = REVISION.parent
REPORT   = (GENETICS.parent.parent / "Platypus - Genomics - General" / "DaRT Data" /
            "DPla24-9425" / "Report-DPla24-9425" /
            "Report_DPla24-9425_1_moreOrders_SNP_mapping_2.csv")
for d in (DERIVED, FIGS):
    d.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------- house style
plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":          8.5,
    "axes.labelsize":     9.0,
    "axes.titlesize":     9.5,
    "axes.titleweight":   "bold",
    "xtick.labelsize":    8.0,
    "ytick.labelsize":    8.0,
    "legend.fontsize":    7.8,
    "legend.frameon":     False,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.linewidth":     0.8,
    "xtick.major.width":  0.8,
    "ytick.major.width":  0.8,
    "lines.linewidth":    1.3,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    "mathtext.default":   "regular",
})

# Okabe-Ito, colourblind-safe
C = {"BORDER": "#0072B2", "SNOWY": "#009E73", "U_MURRAY": "#CC79A7", "KI": "#D55E00"}
M = {"BORDER": "o",       "SNOWY": "s",       "U_MURRAY": "^",       "KI": "D"}
REGION_LABEL = {"BORDER":   "Border Rivers (NSW)",
                "SNOWY":    "Snowy (NSW)",
                "U_MURRAY": "Upper Murray (Vic.)",
                "KI":       "Kangaroo Island (S.A.)"}
GREY = "#4d4d4d"

HO = r"$H_\mathrm{O}$"; HE = r"$H_\mathrm{E}$"
FIS = r"$F_\mathrm{IS}$"; FST = r"$F_\mathrm{ST}$"; AR = r"$A_\mathrm{r}$"

PRETTY = {"KI": "Kangaroo Island", "EUCUMBENE_ABOVE": "Eucumbene above",
          "EUCUMBENE_BELOW": "Eucumbene below", "MITTA_ABOVE": "Mitta above",
          "MITTA_BELOW": "Mitta below", "OVENS": "Ovens",
          "SEVERN_ABOVE": "Severn above", "SEVERN_BELOW": "Severn below",
          "SNOWY": "Snowy", "TENTERFIELD": "Tenterfield", "THREDBO": "Thredbo"}


def save(fig, stem):
    png, pdf = FIGS / f"{stem}.png", FIGS / f"{stem}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"  wrote {png.name} ({png.stat().st_size:,} B) and {pdf.name} ({pdf.stat().st_size:,} B)")


# ------------------------------------------------------- 95% confidence ellipse
# chi-square 0.95 quantile on 2 degrees of freedom, hard-coded so that no scipy
# dependency is introduced (scipy is not installed in the analysis environment).
CHI2_95_DF2 = 5.991464547107979


def conf_ellipse(ax, x, y, colour, lw=1.1, ls="-", alpha=0.95, zorder=3, label=None):
    """Draw the 95% ellipse of a bivariate normal fitted to (x, y).

    The ellipse is obtained from the eigen-decomposition of the sample
    covariance matrix (ddof = 1): the semi-axes are sqrt(chi2_{0.95,2} * lambda_i)
    and the orientation is the angle of the leading eigenvector.  This is the
    standard 95% concentration ellipse — the contour of constant Mahalanobis
    distance that contains 95% of the fitted distribution — not a confidence
    region for the group mean.  It needs at least three points.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    if x.size < 3:
        return None
    cov = np.cov(np.vstack((x, y)), ddof=1)
    if not np.all(np.isfinite(cov)):
        return None
    vals, vecs = np.linalg.eigh(cov)          # ascending
    vals = vals[::-1]; vecs = vecs[:, ::-1]   # leading first
    vals = np.clip(vals, 0.0, None)
    if vals[0] <= 0:
        return None
    w, h = 2.0 * np.sqrt(CHI2_95_DF2 * vals)  # full widths
    ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    e = Ellipse((x.mean(), y.mean()), width=w, height=h, angle=ang,
                facecolor="none", edgecolor=colour, lw=lw, ls=ls, alpha=alpha,
                zorder=zorder, label=label)
    ax.add_patch(e)
    return e


# =============================================================== load out/
print("loading out/ ...")
pop   = pd.read_csv(OUT / "diversity_by_population.csv")
ind   = pd.read_csv(OUT / "diversity_by_individual.csv")
pca   = pd.read_csv(OUT / "pca_scores.csv")
samp  = pd.read_csv(OUT / "samples_used.csv")
fstp  = pd.read_csv(OUT / "fst_pairwise.csv")
fstm  = pd.read_csv(OUT / "fst_matrix.csv", index_col=0)
oldnew= pd.read_csv(OUT / "TableS_old_vs_new.csv")
kinp  = pd.read_csv(OUT / "kinship_pairs.csv")

J   = lambda f: json.loads((OUT / f).read_text())
sfs        = J("sfs.json")
rare       = J("rarefaction.json")
ar         = J("allelic_richness.json")
kin        = J("kinship_summary.json")
polym      = J("polymorphism_matched_n.json")
eigen      = J("pca_eigen.json")
cascade    = J("filter_cascade.json")
lclass     = J("locus_class_recovery.json")
mafsens    = J("maf_sensitivity.json")
ratios     = J("old_vs_new_ratios.json")
locbias    = J("oldnew_locus_bias.json")
jsum       = J("joint_summary.json")

REG = dict(zip(pop.population, pop.region))
REG["KI"] = "KI"
ind["region"] = ind.population.map(REG)
KI_ORDER = ["KI1", "KI2", "KI3", "KI4", "KI5", "KI6", "KI7", "KI8"]
# Sex and age class at capture, from Hawke et al. (2025) Australian Mammalogy 47, AM24042,
# Table A2 (capture year 2021). Field observations, not derived from the genotype data.
KI_SEX_AGE = {
    "KI1": ("Male", "Juvenile"),  "KI2": ("Female", "Juvenile"),
    "KI3": ("Female", "Adult"),   "KI4": ("Male", "Adult"),
    "KI5": ("Male", "Sub-adult"), "KI6": ("Female", "Adult"),
    "KI7": ("Male", "Juvenile"),  "KI8": ("Female", "Juvenile"),
}
VAR = eigen["variance_pct"]


# ============================================ derived: full-n Ovens folded SFS
def ovens_full_spectrum():
    """Folded SFS of the full n = 19 Ovens sample, by minor allele count out of
    38 gene copies.  out/ keeps only its summary statistics, so it is recomputed
    here from the DArT co-processed report on the retained locus set and checked
    against those stored statistics before use."""
    cache = DERIVED / "ovens_full_sfs.json"
    if cache.exists():
        return json.loads(cache.read_text())
    if not REPORT.exists():
        print("  ! DArT report not found; full-n Ovens spectrum unavailable")
        return None
    s = samp[(samp.population == "OVENS") & (samp.included == 1)]
    cols = [0] + [26 + int(c) for c in s.col_index]
    df = pd.read_csv(REPORT, skiprows=6, usecols=cols, low_memory=False)
    keep = set(pd.read_csv(OUT / "loci_final.csv",
                           usecols=["AlleleID", "retained"]).query("retained == 1").AlleleID)
    d = df[df.iloc[:, 0].isin(keep)]
    raw = d.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    # DArT one-row SNP report coding: 0 = hom reference, 1 = hom SNP, 2 = heterozygote
    dose = np.where(raw == 1, 2.0, np.where(raw == 2, 1.0, raw))
    n = dose.shape[1]
    p = dose.sum(axis=1) / (2 * n)
    f = np.minimum(p, 1 - p)
    poly = f > 0
    mac = np.round(f[poly] * 2 * n).astype(int)
    counts = np.bincount(mac, minlength=n + 1)[1:n + 1]
    exp_np, exp_mf = sfs["OVENS_full"]["n_polymorphic"], sfs["OVENS_full"]["mean_folded_maf"]
    assert int(poly.sum()) == exp_np, (poly.sum(), exp_np)
    assert abs(float(f[poly].mean()) - exp_mf) < 1e-9, (f[poly].mean(), exp_mf)
    out = {"n": n, "gene_copies": 2 * n, "n_polymorphic": int(poly.sum()),
           "mean_folded_maf": float(f[poly].mean()),
           "mac_counts": {str(i + 1): int(c) for i, c in enumerate(counts)},
           "mac_proportions": {str(i + 1): float(c / poly.sum()) for i, c in enumerate(counts)},
           "validated_against": {"out/sfs.json OVENS_full n_polymorphic": exp_np,
                                 "out/sfs.json OVENS_full mean_folded_maf": exp_mf}}
    cache.write_text(json.dumps(out, indent=1))
    print(f"  recomputed + validated full-n Ovens spectrum -> {cache.name}")
    return out


# =================================================================== FIGURE 1
def figure1():
    """Figure 1 is *copied*, not regenerated — see the warning printed below.

    Figure 1 is built by <Genetics>/build_fig1_final.py, which needs
      (i)  the pyshp package (`import shapefile`), and
      (ii) the Geoscience Australia coastline/state-boundary shapefile
           `GIS DATA/australia_bndGDA94.shp`, used for BOTH panels.
    Neither is present in this environment (the shapefile is not anywhere under
    the mounted drives; only the KI-local RockyRiver / Platypus_Survey_Sites
    layers are), so the figure cannot be rebuilt here and is copied as-is.

    CONTENT OF THE COPIED IMAGE IS CORRECT (checked 2026-09-14 against
    Hawke et al. 2025, Australian Mammalogy 47, AM24042):
      * panel (b) title      "... capture sites of the 8 genotyped individuals (2021)"
      * panel (b) legend     "Capture site of genotyped individuals\\n(2021; n=8 across 4 sites)"
    Table A2 of that paper lists KI1-KI8 under capture year 2021, and its Results
    state that eight platypuses were captured in May 2021 while March 2022 gave
    six captures of five individuals, one a 2021 recapture.  So the "(2021)" year
    label is right, and build_fig1_final.py's selection of sites on the
    `Platypus21` attribute of Platypus_Survey_Sites.dbf alone (2021 captures:
    8 animals at 4 sites) is the correct selection for the genotyped sample.
    The `Platypus22` attribute (6 further animals at 6 sites, `Sum21_22` = 14
    over both years) covers animals that were not genotyped.  An earlier note
    here claimed the label and the site set were wrong; that claim was mistaken
    and is withdrawn.

    REMAINING COSMETIC ISSUE (cannot be fixed here): the region colours in
    panel (a) - Border Rivers green, Snowy purple, Upper Murray blue - do not
    match the Okabe-Ito scheme used for the same regions in Figure 2 (Border
    Rivers blue, Snowy green, Upper Murray pink), although the Figure 1 caption
    says the regions are grouped as in the ordination.  Redraw on the Figure 2
    palette when the shapefile and pyshp are available.
    """
    print("Figure 1 (copy) ...")
    print("  ! WARNING: Figure 1 is COPIED, not regenerated — pyshp and")
    print("  !          'GIS DATA/australia_bndGDA94.shp' are unavailable here.")
    print("  !          Its content is correct: all eight genotyped animals were")
    print("  !          captured in May 2021 (Hawke et al. 2025, Table A2), so the")
    print("  !          '(2021)' label and the 4-site selection are right as drawn.")
    print("  !          Outstanding cosmetic issue only: panel (a) region colours do")
    print("  !          not match Figure 2's Okabe-Ito scheme. See the docstring.")
    ok = True
    for ext in ("png", "pdf"):
        src = GENETICS / f"Fig1_KI_map_REVISED.{ext}"
        dst = FIGS / f"Figure_1_study_area.{ext}"
        if not src.exists():
            print(f"  ! missing {src}"); ok = False; continue
        shutil.copyfile(src, dst)
        print(f"  wrote {dst.name} ({dst.stat().st_size:,} B)")
    return ok


# =================================================================== FIGURE 2
def panel_tag(ax, tag, x=0.0, y=1.015):
    """Bare panel identifier, outside the plotting area at the top-left corner.

    House style: no titles, no subheadings, no in-panel annotations.  The panel
    identifier is a bare letter in parentheses only, so that the captions can
    refer to '(a)' and '(b)'; everything the figure used to say in words now
    lives in the caption.
    """
    ax.text(x, y, tag, transform=ax.transAxes, fontsize=9.5, fontweight="bold",
            ha="left", va="bottom", clip_on=False)


def figure2():
    print("Figure 2 PCA ...")
    fig = plt.figure(figsize=(7.4, 5.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.19], hspace=0.50, wspace=0.26)

    ki   = pca[pca.population == "KI"]
    main = pca[pca.population != "KI"]

    def scatter(ax, xk, yk):
        for reg in ("BORDER", "SNOWY", "U_MURRAY"):
            g = main[main.region == reg]
            ax.scatter(g[xk], g[yk], s=17, marker=M[reg], facecolor=C[reg],
                       edgecolor="white", linewidth=0.25, alpha=0.85, zorder=2,
                       label=REGION_LABEL[reg])
            # 95% ellipse for the regional group (n = 79, 99, 36)
            conf_ellipse(ax, g[xk], g[yk], C[reg], lw=1.0, ls="-", alpha=0.9, zorder=3)
        ax.scatter(ki[xk], ki[yk], s=52, marker=M["KI"], facecolor=C["KI"],
                   edgecolor="black", linewidth=0.6, zorder=5,
                   label=REGION_LABEL["KI"])
        # 95% ellipse for KI; dashed, because n = 8 makes it imprecise
        conf_ellipse(ax, ki[xk], ki[yk], C["KI"], lw=1.2, ls=(0, (4, 2)),
                     alpha=0.95, zorder=4)
        ax.axhline(0, color="#cccccc", lw=0.6, zorder=0)
        ax.axvline(0, color="#cccccc", lw=0.6, zorder=0)

    # ---- panel (a) PC1 x PC2
    axA = fig.add_subplot(gs[0, 0])
    lo, hi = eigen["overlap"]["PC1"]["mainland_range"]
    axA.axvspan(lo, hi, color="#bbbbbb", alpha=0.20, zorder=0, lw=0)
    scatter(axA, "PC1", "PC2")
    axA.set_xlabel(f"PC1 ({VAR[0]:.2f}% of variance)")
    axA.set_ylabel(f"PC2 ({VAR[1]:.2f}% of variance)")
    panel_tag(axA, "(a)")

    # ---- panel (b) PC2 x PC3
    axB = fig.add_subplot(gs[0, 1])
    scatter(axB, "PC2", "PC3")
    axB.set_xlabel(f"PC2 ({VAR[1]:.2f}% of variance)")
    axB.set_ylabel(f"PC3 ({VAR[2]:.2f}% of variance)")
    panel_tag(axB, "(b)")

    # ---- marginal strips (graphical only; what they show is stated in the caption)
    rng = np.random.default_rng(42)

    def strip(ax, key, parent):
        for reg in ("BORDER", "SNOWY", "U_MURRAY"):
            g = main[main.region == reg]
            ax.scatter(g[key], rng.uniform(-0.30, 0.30, len(g)), s=9, marker=M[reg],
                       facecolor=C[reg], edgecolor="none", alpha=0.65, zorder=2)
        ax.scatter(ki[key], rng.uniform(-0.30, 0.30, len(ki)), s=34, marker=M["KI"],
                   facecolor=C["KI"], edgecolor="black", linewidth=0.5, zorder=5)
        ax.set_ylim(-1.0, 1.0)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlim(parent.get_xlim())
        ax.set_xlabel(key, fontsize=8.0, labelpad=2)
        ax.tick_params(labelsize=7.2)

    axAm = fig.add_subplot(gs[1, 0])
    strip(axAm, "PC1", axA)
    axAm.axvspan(lo, hi, color="#bbbbbb", alpha=0.20, zorder=0, lw=0)

    axBm = fig.add_subplot(gs[1, 1])
    strip(axBm, "PC2", axB)

    handles = [Line2D([0], [0], marker=M[r], color="none", markerfacecolor=C[r],
                      markeredgecolor="white", markersize=6, label=REGION_LABEL[r])
               for r in ("BORDER", "SNOWY", "U_MURRAY")]
    handles.append(Line2D([0], [0], marker="D", color="none", markerfacecolor=C["KI"],
                          markeredgecolor="black", markersize=7.5,
                          label=REGION_LABEL["KI"]))
    fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.055),
               handletextpad=0.5, columnspacing=1.8)
    save(fig, "Figure_2_PCA")


# =================================================================== FIGURE 3
def figure3():
    print("Figure 3 individual diversity ...")
    order = ["KI"] + list(pop[pop.population != "KI"].sort_values("Ho", ascending=False).population)
    xpos = {p: i for i, p in enumerate(order)}
    rng = np.random.default_rng(7)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.6, 4.6))
    fig.subplots_adjust(wspace=0.30)

    ki_max_ho = float(ind.loc[ind.population == "KI", "Ho_ind"].max())

    # ---- (a) population H_O with individuals overlaid
    for p in order:
        x = xpos[p]; reg = REG[p]
        g = ind[ind.population == p]
        axA.scatter(x + rng.uniform(-0.22, 0.22, len(g)), g.Ho_ind, s=11,
                    marker=M[reg], facecolor=C[reg], edgecolor="none", alpha=0.55, zorder=2)
        r = pop[pop.population == p].iloc[0]
        axA.errorbar(x, r.Ho, yerr=r.Ho_SE, fmt="_", color="black", markersize=17,
                     markeredgewidth=1.7, elinewidth=1.1, capsize=2.6, zorder=6)
        axA.plot([x - 0.30, x + 0.30], [r.He_unbiased] * 2, color=GREY, lw=1.1,
                 ls=(0, (3, 2)), zorder=5)
    axA.axhline(ki_max_ho, color=C["KI"], lw=0.9, ls=":", zorder=1)
    axA.set_xticks(range(len(order)))
    axA.set_xticklabels([PRETTY[p] for p in order], rotation=48, ha="right")
    axA.set_ylabel("Observed heterozygosity, " + HO)
    axA.set_xlim(-0.7, len(order) - 0.3)
    panel_tag(axA, "(a)")

    # ---- (b) individual F in the mainland reference frame
    mmax = jsum["individuals"]["mainland_F_uni_mainland_frame_max"]
    mmean = jsum["individuals"]["mainland_F_uni_mainland_frame_mean"]
    for p in order:
        x = xpos[p]; reg = REG[p]
        g = ind[ind.population == p]
        axB.scatter(x + rng.uniform(-0.22, 0.22, len(g)), g.F_uni_mainland, s=11,
                    marker=M[reg], facecolor=C[reg], edgecolor="none",
                    alpha=0.85 if p == "KI" else 0.55, zorder=3 if p == "KI" else 2)
    axB.axhline(mmax, color=GREY, lw=0.9, ls="--", zorder=1)
    axB.axhline(mmean, color=GREY, lw=0.7, ls=":", zorder=1)
    axB.axhline(0, color="#cccccc", lw=0.7, zorder=0)
    axB.set_xticks(range(len(order)))
    axB.set_xticklabels([PRETTY[p] for p in order], rotation=48, ha="right")
    axB.set_ylabel(r"Individual inbreeding $F_\mathrm{UNI}$ (mainland frame)")
    axB.set_xlim(-0.7, len(order) - 0.3)
    panel_tag(axB, "(b)")

    handles = [
        Line2D([0], [0], marker="_", color="black", lw=0, markersize=13,
               markeredgewidth=1.7, label="population " + HO + " ± SE across loci"),
        Line2D([0], [0], color=GREY, lw=1.1, ls=(0, (3, 2)), label="population " + HE),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#999999",
               markersize=4.5, label="individual value"),
        Line2D([0], [0], color=C["KI"], lw=0.9, ls=":", label="highest KI individual, (a)"),
        Line2D([0], [0], color=GREY, lw=0.9, ls="--", label="highest mainland individual, (b)"),
        Line2D([0], [0], color=GREY, lw=0.7, ls=":", label="mainland mean, (b)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.20),
               handletextpad=0.7, columnspacing=2.0, labelspacing=0.45)
    save(fig, "Figure_3_diversity_individual")


# =================================================================== FIGURE 4
def figure4(ovens_full):
    print("Figure 4 folded SFS at matched n ...")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.5),
                                   gridspec_kw=dict(width_ratios=[1.45, 1.0], wspace=0.62))

    mac = [str(i) for i in range(1, 9)]
    x = np.arange(1, 9) / 16.0
    ki_p = np.array([sfs["KI"]["mac_proportions"][k] for k in mac])
    env_pops = ["OVENS", "SNOWY", "THREDBO", "TENTERFIELD"]
    env = np.array([[sfs[f"{p}_rarefied_n8"]["mac_proportions_mean"][k] for k in mac]
                    for p in env_pops])

    axA.fill_between(x, env.min(0), env.max(0), color=C["SNOWY"], alpha=0.22, lw=0,
                     zorder=1, label="mainland envelope, rarefied to n = 8")
    axA.plot(x, env.mean(0), color=C["SNOWY"], lw=1.5, marker="s", ms=4.2, zorder=3,
             label="mainland mean, n = 8")
    axA.plot(x, ki_p, color=C["KI"], lw=1.8, marker="D", ms=5.0, zorder=4,
             label="Kangaroo Island, n = 8")
    handles_extra = []
    if ovens_full:
        nfull = ovens_full["n"]
        xf = np.arange(1, nfull + 1) / (2 * nfull)
        yf = np.array([ovens_full["mac_proportions"][str(i)] for i in range(1, nfull + 1)])
        axA.plot(xf, yf, color=GREY, lw=1.0, ls="--", marker=".", ms=3.0, zorder=2,
                 label="Ovens, full n = 19 (unmatched)")
        axA.axvline(1 / 16, color="#999999", lw=0.7, ls=":", zorder=0)
        handles_extra.append(Line2D([0], [0], color="#999999", lw=0.7, ls=":",
                                    label="smallest folded MAF at n = 8 (1/16)"))
    axA.set_xlabel("Folded minor allele frequency")
    axA.set_ylabel("Proportion of polymorphic loci")
    panel_tag(axA, "(a)")

    # ---- (b) decomposition of the published KI-Ovens gap
    vals = [("Ovens, full n = 19\n(unmatched)", sfs["OVENS_full"]["mean_folded_maf"], GREY),
            ("mainland, rarefied\nto n = 8", float(np.mean([sfs[f"{p}_rarefied_n8"]["mean_folded_maf"]
                                                            for p in env_pops])), C["SNOWY"]),
            ("Kangaroo Island,\nn = 8", sfs["KI"]["mean_folded_maf"], C["KI"])]
    ypos = np.arange(len(vals))[::-1]
    for y, (lab, v, col) in zip(ypos, vals):
        axB.barh(y, v, height=0.52, color=col, alpha=0.9, zorder=2)
    axB.set_yticks(ypos)
    axB.set_yticklabels([v[0] for v in vals], fontsize=7.4, linespacing=1.25)
    axB.set_xlim(0, 0.30)
    axB.set_ylim(-0.6, len(vals) - 0.4)
    axB.set_xlabel("Mean folded minor allele frequency")
    panel_tag(axB, "(b)")

    h, _ = axA.get_legend_handles_labels()
    fig.legend(handles=h + handles_extra, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.16), handletextpad=0.7, columnspacing=1.8,
               labelspacing=0.45)
    save(fig, "Figure_4_SFS_matched_n")


# ================================================================== FIGURE S1
def figureS1():
    print("Figure S1 rarefaction ...")
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 3.0))
    fig.subplots_adjust(wspace=0.40)
    axA, axB, axC = axes
    levels = [5, 8, 10, 15, 20]
    kiobs = rare["KI_observed"]

    for p, d in rare["by_pop"].items():
        reg = REG[p]
        ns = [n for n in levels if f"n{n}" in d] + [d["n_full"]]
        he = [d[f"n{n}"]["He_mean"] for n in levels if f"n{n}" in d] + [d["full"]["He"]]
        sd = [d[f"n{n}"]["He_sd"] for n in levels if f"n{n}" in d] + [0.0]
        axA.errorbar(ns, he, yerr=sd, color=C[reg], marker=M[reg], ms=3.4, lw=1.0,
                     elinewidth=0.8, capsize=2, alpha=0.9)
        npl = [d[f"n{n}"]["n_poly_mean"] for n in levels if f"n{n}" in d] + [d["full"]["n_poly"]]
        nsd = [d[f"n{n}"]["n_poly_sd"] for n in levels if f"n{n}" in d] + [0.0]
        axB.errorbar(ns, npl, yerr=nsd, color=C[reg], marker=M[reg], ms=3.4, lw=1.0,
                     elinewidth=0.8, capsize=2, alpha=0.9)

    for ax, val in ((axA, kiobs["He"]), (axB, kiobs["n_poly"])):
        ax.axhline(val, color=C["KI"], ls="--", lw=1.1, zorder=1)
        ax.scatter([8], [val], marker="D", s=46, color=C["KI"], edgecolor="black",
                   linewidth=0.6, zorder=6)
    axA.set_xlabel("Subsample size (individuals)"); axA.set_ylabel("Expected heterozygosity, " + HE)
    axB.set_xlabel("Subsample size (individuals)"); axB.set_ylabel("Polymorphic loci")
    panel_tag(axA, "(a)"); panel_tag(axB, "(b)")

    gs_ = [2, 4, 6, 8, 10, 12, 14, 16]
    for p, cur in ar["curve"].items():
        reg = REG[p]
        y = [cur[str(g)] for g in gs_]
        if p == "KI":
            axC.plot(gs_, y, color=C["KI"], lw=1.9, marker="D", ms=4.0, zorder=6)
        else:
            axC.plot(gs_, y, color=C[reg], lw=0.9, marker=M[reg], ms=2.8, alpha=0.75)
    axC.set_xlabel("Gene copies sampled, g")
    axC.set_ylabel("Rarefied allelic richness, " + AR)
    panel_tag(axC, "(c)")

    handles = [Line2D([0], [0], color=C[r], marker=M[r], ms=4, lw=1.1, label=REGION_LABEL[r])
               for r in ("BORDER", "SNOWY", "U_MURRAY")]
    handles.append(Line2D([0], [0], color=C["KI"], marker="D", ms=5, lw=1.6,
                          label=REGION_LABEL["KI"]))
    handles.append(Line2D([0], [0], color=C["KI"], ls="--", lw=1.1,
                          label="observed KI value, (a) and (b)"))
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.16),
               handletextpad=0.5, columnspacing=1.6, labelspacing=0.45)
    save(fig, "Figure_S1_rarefaction")


# ================================================================== FIGURE S2
def figureS2():
    print("Figure S2 kinship ...")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.3),
                                   gridspec_kw=dict(width_ratios=[1.25, 1.0], wspace=0.30))
    six = kin["matched_n8_null"]["pops"]
    full_null = kinp[(kinp.within_pop == 1) & (kinp.population_i.isin(six))]["relatedness_scale"].to_numpy()
    ki_vals = np.array(kin["KI"]["values"], float)
    m8 = kin["matched_n8_null"]["kinship"]

    axA.hist(full_null, bins=60, color="#b8b8b8", edgecolor="none", zorder=1,
             label="mainland within-population pairs, full sample")
    axA.axvspan(m8["p05"], m8["p95"], color=C["SNOWY"], alpha=0.20, lw=0, zorder=2,
                label="matched n = 8 null, 5th–95th percentile")
    extra = []
    for v, ls, lab in ((m8["p95"], "--", "95th percentile"), (m8["p99"], "-.", "99th percentile"),
                       (m8["max"], ":", "maximum")):
        axA.axvline(v, color=C["SNOWY"], lw=1.0, ls=ls, zorder=3)
        extra.append(Line2D([0], [0], color=C["SNOWY"], lw=1.0, ls=ls,
                            label=f"matched n = 8 null, {lab}"))
    yr = axA.get_ylim()[1]
    axA.set_ylim(0, yr * 1.30)
    axA.scatter(ki_vals, np.full_like(ki_vals, yr * 1.13), marker="D", s=26,
                color=C["KI"], edgecolor="black", linewidth=0.4, zorder=6,
                label="the 28 KI pairs")
    axA.set_xlabel("Pairwise relatedness (Yang et al. 2010)")
    axA.set_ylabel("Number of pairs")
    panel_tag(axA, "(a)")

    o = np.argsort(ki_vals)
    axB.axhspan(m8["p05"], m8["p95"], color=C["SNOWY"], alpha=0.20, lw=0, zorder=1)
    axB.axhline(m8["p95"], color=C["SNOWY"], lw=1.0, ls="--", zorder=2)
    axB.axhline(m8["max"], color=C["SNOWY"], lw=1.0, ls=":", zorder=2)
    axB.axhline(0.125, color=GREY, lw=0.9, ls="-.", zorder=2)
    axB.scatter(np.arange(28), ki_vals[o], marker="D", s=24, color=C["KI"],
                edgecolor="black", linewidth=0.4, zorder=5)
    axB.set_xlabel("KI pair, ranked")
    axB.set_ylabel("Pairwise relatedness")
    axB.set_xlim(-1, 29)
    axB.set_ylim(-0.52, 0.40)
    panel_tag(axB, "(b)")

    h, _ = axA.get_legend_handles_labels()
    h = h + extra + [Line2D([0], [0], color=GREY, lw=0.9, ls="-.",
                            label="nominal second-order threshold (0.125)")]
    fig.legend(handles=h, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.28),
               handletextpad=0.7, columnspacing=2.0, labelspacing=0.45)
    save(fig, "Figure_S2_kinship")


# ================================================================== FIGURE S3
def figureS3(ovens_full):
    print("Figure S3 legacy vs matched SFS ...")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True,
                                   gridspec_kw=dict(wspace=0.13))
    mac = [str(i) for i in range(1, 9)]
    x8 = np.arange(1, 9) / 16.0
    ki_p = np.array([sfs["KI"]["mac_proportions"][k] for k in mac])

    axA.plot(x8, ki_p, color=C["KI"], lw=1.8, marker="D", ms=5.0, zorder=4,
             label="Kangaroo Island, n = 8")
    if ovens_full:
        nf = ovens_full["n"]
        xf = np.arange(1, nf + 1) / (2 * nf)
        yf = [ovens_full["mac_proportions"][str(i)] for i in range(1, nf + 1)]
        axA.plot(xf, yf, color=GREY, lw=1.3, marker="o", ms=3.2, zorder=3,
                 label="Ovens, full n = 19")
    axA.axvline(1 / 16, color="#999999", lw=0.7, ls=":", zorder=0)
    axA.set_xlabel("Folded minor allele frequency")
    axA.set_ylabel("Proportion of polymorphic loci")
    axA.set_ylim(-0.015, 0.40)
    panel_tag(axA, "(a)")

    envm = np.array([sfs["OVENS_rarefied_n8"]["mac_proportions_mean"][k] for k in mac])
    envlo = np.array([sfs["OVENS_rarefied_n8"]["mac_proportions_lo"][k] for k in mac])
    envhi = np.array([sfs["OVENS_rarefied_n8"]["mac_proportions_hi"][k] for k in mac])
    axB.fill_between(x8, envlo, envhi, color=C["SNOWY"], alpha=0.25, lw=0, zorder=1,
                     label="Ovens n = 8, 95% range of 200 subsamples")
    axB.plot(x8, envm, color=C["SNOWY"], lw=1.5, marker="s", ms=4.2, zorder=3,
             label="Ovens, rarefied to n = 8")
    axB.plot(x8, ki_p, color=C["KI"], lw=1.8, marker="D", ms=5.0, zorder=4,
             label="Kangaroo Island, n = 8")
    axB.set_xlabel("Folded minor allele frequency")
    panel_tag(axB, "(b)")

    hA, lA = axA.get_legend_handles_labels()
    hB, lB = axB.get_legend_handles_labels()
    handles, labels = [], []
    for hh, ll in zip(hA + hB, lA + lB):
        if ll not in labels:
            handles.append(hh); labels.append(ll)
    handles.append(Line2D([0], [0], color="#999999", lw=0.7, ls=":",
                          label="smallest folded MAF at n = 8 (1/16)"))
    fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.24),
               handletextpad=0.7, columnspacing=2.0, labelspacing=0.45)
    save(fig, "Figure_S3_SFS_legacy_vs_matched")


# ===================================================================== TABLES
def f4(v):  return "—" if pd.isna(v) else f"{v:.4f}"
def f3(v):  return "—" if pd.isna(v) else f"{v:.3f}"
def s4(v):  return "—" if pd.isna(v) else f"{v:+.4f}"
def s3(v):  return "—" if pd.isna(v) else f"{v:+.3f}"
def pc1(v): return "—" if pd.isna(v) else f"{v:.1f}%"
def ithou(v): return "—" if pd.isna(v) else f"{int(round(v)):,}"


def tables():
    print("TABLES.md ...")
    L = []
    A = L.append

    A("# Tables for the revised manuscript")
    A("")
    A("*Drifting alone: genome-wide diversity in the isolated, introduced Kangaroo Island "
      "platypus (*Ornithorhynchus anatinus*)* — Bino, Hawke, Baring & Gongora.")
    A("")
    A("All values are from the co-processed re-analysis (`work/out/`, generated "
      "2026-09-12 by `ki_coprocessed_joint.py`, seed 42). Numbering follows "
      "REVISION\\_PLAN.md Part C2. Every table uses the same filtered locus set: "
      "**4,002 autosomal SNPs genotyped in all 222 individuals (214 mainland + 8 KI) "
      "with zero missing data**.")
    A("")
    A("---")
    A("")

    # ---------------- Table 1
    o = ["KI"] + [p for p in pop.population if p != "KI"]
    d = pop.set_index("population").loc[o].reset_index()
    A("## Table 1. Population-level diversity")
    A("")
    A("| Population | Region | *n* | H_O ± SE | H_E ± SE | F_IS | jackknife SE | 95% bootstrap CI | "
      "Polymorphic loci | Polymorphic (%) | A_r (g = 16) |")
    A("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for _, r in d.iterrows():
        name = "**Kangaroo Island**" if r.population == "KI" else PRETTY[r.population]
        A(f"| {name} | {REGION_LABEL[r.region].split(' (')[0]} | {int(r.n)} | "
          f"{f4(r.Ho)} ± {f4(r.Ho_SE)} | {f4(r.He_unbiased)} ± {f4(r.He_SE)} | "
          f"{s4(r.Fis_WC)} | {f4(r.Fis_jack_SE)} | [{s4(r.Fis_boot_lo)}, {s4(r.Fis_boot_hi)}] | "
          f"{int(r.n_polymorphic):,} | {pc1(100*r.prop_polymorphic)} | {f4(r.Ar_g16)} |")
    mn = pop[pop.population != "KI"]
    A(f"| *Mainland mean (10 populations)* | — | 214 | "
      f"{f4(mn.Ho.mean())} ± {f4(mn.Ho.std(ddof=1)/np.sqrt(len(mn)))} | "
      f"{f4(mn.He_unbiased.mean())} ± {f4(mn.He_unbiased.std(ddof=1)/np.sqrt(len(mn)))} | "
      f"{s4(mn.Fis_WC.mean())} | — | — | "
      f"{ithou(mn.n_polymorphic.mean())} | {pc1(100*mn.prop_polymorphic.mean())} | "
      f"{f4(mn.Ar_g16.mean())} |")
    A("")
    A("**Table 1.** Genome-wide diversity of the Kangaroo Island (KI) platypus and ten "
      "mainland Australian platypus populations, computed on a single co-processed DArTseq "
      "SNP panel of 4,002 autosomal loci scored in all 222 individuals. "
      "H_O, observed heterozygosity; H_E, unbiased expected heterozygosity (Nei 1987); "
      "F_IS, Weir & Cockerham (1984) inbreeding coefficient as a ratio of sums over loci; "
      "A_r, rarefied allelic richness at g = 16 gene copies (El Mousadik & Petit 1996). "
      "**± values for H_O and H_E are standard errors across loci** (the standard deviation "
      "of the per-locus estimates divided by the square root of the number of loci). "
      "**F_IS carries a delete-one-locus jackknife standard error and a 1,000-replicate "
      "bootstrap-over-loci 95% confidence interval**, not a standard error across loci. "
      "**On the mainland-mean row the ± value is the standard error among the ten population "
      "estimates, not across loci.** A_r cannot be standardised to g = 16 for the two "
      "populations with n = 4 (Eucumbene above, Mitta below), and the mainland mean for A_r "
      "is therefore taken over the eight populations that can be standardised. "
      "Absolute heterozygosities are not comparable with the published Mijangos et al. (2022) "
      "values because the ascertainment panel differs (see Table S3); all KI-versus-mainland "
      "statements are made in relative terms.")
    A("")
    A(f"KI H_E is {ratios['KI_pct_of_mainland_mean_He_new']:.1f}% of the mainland mean "
      f"({f4(mn.He_unbiased.mean())}) and {ratios['KI_pct_of_OVENS_He_new']:.1f}% of the Ovens "
      f"value; KI is the lowest of the eleven populations for H_E, H_O and A_r.")
    A("")
    A("---")
    A("")

    # ---------------- Table 2
    kid = ind[ind.population == "KI"].set_index("sample").loc[KI_ORDER].reset_index()
    panel_mean_ho = float(ind.Ho_ind.mean())
    kid["sMLH"] = kid.Ho_ind / panel_mean_ho
    kid["sex"] = [KI_SEX_AGE[s][0] for s in kid["sample"]]
    kid["age_class"] = [KI_SEX_AGE[s][1] for s in kid["sample"]]
    kid.to_csv(DERIVED / "table2_KI_individuals.csv", index=False)
    A("## Table 2. Individual-level heterozygosity and inbreeding, KI1–KI8")
    A("")
    A("| Individual | Sex | Age class | Loci scored | H_O | F_UNI (KI frame) | 95% CI | "
      "F_UNI (mainland frame) | 95% CI |")
    A("|---|---|---|---:|---:|---:|---:|---:|---:|")
    for _, r in kid.iterrows():
        A(f"| {r['sample']} | {KI_SEX_AGE[r['sample']][0]} | "
          f"{KI_SEX_AGE[r['sample']][1]} | {int(r.n_loci_called):,} | {f4(r.Ho_ind)} | "
          f"{s3(r.F_uni_ownpop)} | [{s3(r.F_uni_ownpop_lo95)}, {s3(r.F_uni_ownpop_hi95)}] | "
          f"{s3(r.F_uni_mainland)} | [{s3(r.F_uni_mainland_lo95)}, {s3(r.F_uni_mainland_hi95)}] |")
    A(f"| **KI mean (n = 8)** | — | — | 4,002 | **{f4(kid.Ho_ind.mean())}** | "
      f"{s3(kid.F_uni_ownpop.mean())} | — | "
      f"**{s3(kid.F_uni_mainland.mean())}** | — |")
    A(f"| *KI range* | — | — | — | {f4(kid.Ho_ind.min())}–{f4(kid.Ho_ind.max())} | "
      f"{s3(kid.F_uni_ownpop.min())} to {s3(kid.F_uni_ownpop.max())} | — | "
      f"{s3(kid.F_uni_mainland.min())} to {s3(kid.F_uni_mainland.max())} | — |")
    mi = ind[ind.population != "KI"]
    A(f"| *Mainland, 214 individuals: mean* | — | — | 4,002 | {f4(mi.Ho_ind.mean())} | "
      f"— | — | "
      f"{s3(mi.F_uni_mainland.mean())} | — |")
    A(f"| *Mainland, 214 individuals: range* | — | — | — | {f4(mi.Ho_ind.min())}–{f4(mi.Ho_ind.max())} | "
      f"— | — | "
      f"{s3(mi.F_uni_mainland.min())} to {s3(mi.F_uni_mainland.max())} | — |")
    A("")
    A("**Table 2.** Per-individual metrics for the eight Kangaroo Island platypuses, with "
      "mainland summary rows, on the same 4,002-locus filtered set as Table 1. "
      "H_O, the proportion of scored loci at which the individual is heterozygous. "
      "F_UNI, the correlation of uniting gametes (Yang et al. 2010, F-hat-3), computed in two "
      "reference frames: against KI allele frequencies (1,105 loci polymorphic within KI) and "
      "against the pooled 214-sample mainland allele frequencies (3,414 loci polymorphic in the "
      "mainland panel). **± / bracketed values are 1,000-replicate bootstrap-over-loci 95% "
      "confidence intervals**, not standard errors. Within-sample allele frequencies force "
      "own-frame F towards zero by construction, so KI-frame values describe variation among "
      "the eight animals rather than their absolute level; the mainland frame places all 222 "
      "animals on a common scale. In that frame the KI mean is "
      f"{jsum['individuals']['KI_F_uni_mainland_mean']:+.3f} against a mainland mean of "
      f"{jsum['individuals']['mainland_F_uni_mainland_frame_mean']:+.3f} and a mainland maximum of "
      f"{jsum['individuals']['mainland_F_uni_mainland_frame_max']:+.3f}, which "
      f"{jsum['individuals']['n_KI_above_mainland_max_F_uni']} of the 8 KI animals exceed. "
      "Inbreeding from runs of homozygosity (F_ROH) was not attempted: the marker density of a "
      "reduced-representation SNP panel cannot resolve ROH tracts reliably. "
      "Sex and age class at capture are field observations reported in Table A2 of Hawke et "
      "al. (2025) and are not derived from the genotype data. Capture site is not a column of "
      "this table: all eight animals were captured in the Rocky River, Flinders Chase "
      "National Park, in May 2021, and the individual-to-site assignment is not carried in "
      "the analysis metadata. Standardised multilocus heterozygosity is not reported: with every individual "
      "scored at the same 4,002 loci and no missing data, sMLH is each animal's H_O rescaled "
      f"by the 222-individual panel mean ({f4(panel_mean_ho)}) and adds nothing to the H_O "
      "column. Capture site and sMLH are absent from the manuscript's Table 2 and are absent "
      "here; sex and age class are present in both.")
    A("")
    A("---")
    A("")

    # ---------------- Table 3
    kf = fstp[(fstp.pop_i == "KI") | (fstp.pop_j == "KI")].copy()
    kf["other"] = np.where(kf.pop_i == "KI", kf.pop_j, kf.pop_i)
    kf = kf.sort_values("Fst_WC")
    A("## Table 3. Pairwise differentiation between KI and each mainland population")
    A("")
    A("| Comparison | *n* (mainland population) | F_ST (Weir & Cockerham θ) | 95% bootstrap CI | G_ST (Nei) | G″_ST (Hedrick) |")
    A("|---|---:|---:|---:|---:|---:|")
    for _, r in kf.iterrows():
        n_other = int(r.n_j if r.pop_i == "KI" else r.n_i)
        A(f"| KI vs {PRETTY[r.other]} | {n_other} | **{f3(r.Fst_WC)}** | "
          f"[{f3(r.Fst_lo95)}, {f3(r.Fst_hi95)}] | {f3(r.Gst_Nei)} | {f3(r.GstPP_Hedrick)} |")
    A(f"| *Mean, KI vs mainland* | — | {f3(kf.Fst_WC.mean())} | — | — | — |")
    mm = fstp[(fstp.pop_i != "KI") & (fstp.pop_j != "KI")]
    A(f"| *Mean, mainland vs mainland* | — | {f3(mm.Fst_WC.mean())} | "
      f"range {f3(mm.Fst_WC.min())}–{f3(mm.Fst_WC.max())} | — | — |")
    A("")
    A("**Table 3.** Pairwise genetic differentiation of the Kangaroo Island platypus population "
      "from each of the ten mainland populations, on the 4,002-locus co-processed set. "
      "F_ST is the Weir & Cockerham (1984) θ computed as a ratio of sums over loci. "
      "**Bracketed values are 1,000-replicate bootstrap-over-loci 95% confidence intervals**; "
      "no ± standard errors are reported for this table. G_ST (Nei 1973) and the "
      "standardised G″_ST (Hedrick 2005; Meirmans & Hedrick 2011) are given for comparability "
      "with studies reporting those statistics. Rows are ordered by increasing F_ST: "
      f"KI is least differentiated from the Ovens ({f3(kf.Fst_WC.min())}), the closest "
      "available contemporary analogue of the documented Victorian (Healesville/upper Yarra) "
      f"founder source, and most differentiated from the Severn below ({f3(kf.Fst_WC.max())}). "
      "The full 11 × 11 matrix is given in Table S6.")
    A("")
    A("---")
    A("")

    # ---------------- Table 4
    A("## Table 4. Diversity at matched sample size (rarefaction to n = 8)")
    A("")
    A("| Population | *n* (full) | H_E full | H_E at n = 8 | % of full | Polymorphic loci, full | "
      "Polymorphic loci at n = 8 | A_r (g = 16) at n = 8 | KI as % of this (polymorphic loci) |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    ki_poly = polym["KI"]["n_polymorphic"]
    for p, d in rare["by_pop"].items():
        n8 = d["n8"]
        A(f"| {PRETTY[p]} | {d['n_full']} | {f4(d['full']['He'])} | "
          f"{f4(n8['He_mean'])} ± {f4(n8['He_sd'])} | {pc1(n8['He_pct_of_full'])} | "
          f"{int(d['full']['n_poly']):,} | {ithou(n8['n_poly_mean'])} ± {n8['n_poly_sd']:.0f} | "
          f"{f4(n8['Ar_g16_mean'])} | "
          f"{pc1(polym['mainland_at_n8'][p]['KI_pct_of_this'])} |")
    k = rare["KI_observed"]
    A(f"| **Kangaroo Island (observed)** | **8** | — | **{f4(k['He'])}** | — | — | "
      f"**{int(k['n_poly']):,}** | **{f4(k['Ar_g16'])}** | — |")
    A(f"| *Mainland mean at n = 8* | — | — | — | — | — | "
      f"{ithou(polym['mainland_mean_n_poly_at_n8'])} | "
      f"{f4(np.mean([rare['by_pop'][p]['n8']['Ar_g16_mean'] for p in rare['by_pop']]))} | "
      f"**{polym['KI_pct_of_mainland_mean_at_n8']:.1f}%** |")
    A("")
    A("**Table 4.** Sensitivity of the diversity statistics to sample size. Each mainland "
      "population with n ≥ 19 was subsampled without replacement to n = 8, the Kangaroo Island "
      "sample size, 200 times (seed 42). **± values are standard deviations across the 200 "
      "subsamples**, not standard errors across loci. H_E is essentially unbiased by sample "
      "size — every population recovers within about 1% of its full-sample value at n = 8 — so "
      "the KI H_E deficit is not a sample-size artefact. Counts of polymorphic loci and "
      "allelic richness are strongly sample-size dependent and are compared only at matched n. "
      f"Against the mean of the six matched mainland samples "
      f"({ithou(polym['mainland_mean_n_poly_at_n8'])} loci) KI retains "
      f"{polym['KI_pct_of_mainland_mean_at_n8']:.1f}% of the polymorphic loci, a "
      f"sample-size-controlled deficit of {polym['KI_deficit_vs_mainland_mean_at_n8_pct']:.1f}%.")
    A("")
    A("---")
    A("")
    A("# Supplementary tables")
    A("")

    # ---------------- Table S1
    A("## Table S1. Filter cascade")
    A("")
    A("| Step | Filter | Loci passing (applied alone) | Loci remaining (cumulative) | Removed at this step | Removed (%) |")
    A("|---|---|---:|---:|---:|---:|")
    prev = None
    for s in cascade["steps"]:
        cum = s["cumulative"]
        rem = "—" if prev is None else f"{prev - cum:,}"
        pct = "—" if prev is None else f"{100*(prev-cum)/prev:.1f}%"
        sa = s["standalone"]
        if isinstance(sa, int):
            sa_s = f"{sa:,}"
        elif s["name"].startswith("HWE"):
            sa_s = f"{prev - cum} removed"
        else:
            sa_s = "one per CloneID"
        A(f"| {cascade['steps'].index(s)} | {s['name']} | {sa_s} | **{cum:,}** | {rem} | {pct} |")
        prev = cum
    A("")
    A("**Call-rate sensitivity**")
    A("")
    A("| Minimum call rate over the 222 samples | Final loci | Polymorphic |")
    A("|---|---:|---:|")
    for cr, v in cascade["callrate_sensitivity"].items():
        A(f"| {cr} | {v['final_loci']:,} | {v['polymorphic']:,} |")
    A("")
    A("**Table S1.** Locus filter cascade from the 22,054-locus DArT co-processed report "
      "(both sequencing orders called in a single pipeline run) to the final analysis set. "
      "The \"applied alone\" column gives the number of loci that pass each filter applied to "
      "the full input panel in isolation; the \"cumulative\" column gives the number surviving "
      "that filter and all preceding ones, which is the set carried forward. "
      "Sex-linked loci were removed using a scaffold→chromosome crosswalk (Table S9). "
      f"The final set is **{cascade['final_loci']:,} SNPs × 222 individuals with zero missing "
      f"genotypes** ({cascade['final_loci_polymorphic']:,} polymorphic across the joint panel, "
      f"{cascade['final_loci_monomorphic_in_222']:,} monomorphic), "
      f"{ratios['locus_gain_factor']:.2f}× the 1,296 loci of the submitted analysis. "
      "Hardy–Weinberg filtering was applied per population among the six populations with "
      "n ≥ 19, removing a locus only where it departed in ≥ 2 populations under a Bonferroni "
      f"correction (α = 0.05/(n_loci × 6)); this removed {cascade['hwe']['removed']} loci. "
      "The joint-sample Bonferroni form quoted in the submitted Methods would instead remove "
      f"{cascade['hwe']['sensitivity_joint_bonferroni_removed']} loci and give "
      f"{cascade['hwe']['sensitivity_joint_final_loci']:,} loci; results are unchanged in "
      "substance. The one-SNP-per-CloneID rule retains the SNP with the highest call rate, "
      "then the highest AvgPIC, then the lowest row index. No ± values appear in this table.")
    A("")
    A("---")
    A("")

    # ---------------- Table S2
    A("## Table S2. Locus classes recovered by co-processing")
    A("")
    ip, ff, lost = (lclass["input_panel_22054"], lclass["final_filtered_set"],
                    lclass["final_set_lost_by_old_cloneid_rule"])
    pctd = lclass["pct_of_class_that_the_old_rule_would_have_dropped"]
    rows = [
        ("Polymorphic within DPla19-4171 (Mijangos order)", "polymorphic_in_DPla19-4171_order"),
        ("Polymorphic within DPla24-9425 (2024 order)", "polymorphic_in_DPla24-9425_order"),
        ("Monomorphic in DPla19-4171, polymorphic in DPla24-9425", "monomorphic_in_4171_but_polymorphic_in_9425"),
        ("Monomorphic in DPla24-9425, polymorphic in DPla19-4171", "monomorphic_in_9425_but_polymorphic_in_4171"),
        ("Monomorphic in both orders", "monomorphic_in_both_orders"),
        ("Monomorphic in the mainland 214, polymorphic in KI", "monomorphic_in_Mijangos214_but_polymorphic_in_KI8"),
        ("Polymorphic in the mainland 214, fixed in KI", "polymorphic_in_Mijangos214_but_monomorphic_in_KI8"),
        ("Fixed allelic difference, KI vs mainland", "fixed_allelic_difference_KI_vs_mainland"),
    ]
    A("| Locus class | Input panel (22,054) | Final filtered set (4,002) | Of the final set, dropped by the CloneID rule | Dropped (%) | Recoverable by the CloneID merge? |")
    A("|---|---:|---:|---:|---:|---|")
    for lab, k in rows:
        dr = lost[k]; tot = ff[k]
        p = 100 * dr / tot if tot else float("nan")
        rec = "**No**" if p >= 90 else ("Partly" if p >= 10 else "Yes")
        A(f"| {lab} | {ip[k]:,} | {tot:,} | {dr:,} | {p:.1f}% | {rec} |")
    A("")
    A("**Table S2.** Locus classes present in the co-processed matrix and the extent to which "
      "the CloneID-matching merge used in the submitted manuscript would have discarded them. "
      "Because both DArT orders were called in a single pipeline run, loci that are monomorphic "
      "within one order, and loci that are fixed allelic differences between KI and the "
      "mainland, are present in the matrix and are retained. The CloneID rule — keep a CloneID "
      "only if it is present in both reports and the two SNP call strings match — discarded "
      f"{locbias['lost_by_old']:,} of the {locbias['final_loci']:,} final loci "
      f"({locbias['pct_lost_by_old']:.1f}%), and discarded almost the whole of exactly the "
      "classes that carry the KI signal: 97.8% of loci monomorphic in the Mijangos order but "
      "polymorphic in the 2024 order, 98.2% of loci polymorphic only in KI, and 100% of the "
      "fixed allelic differences between KI and the mainland. "
      "The \"recoverable\" column reads No where the rule would have dropped ≥ 90% of the class, "
      "Partly where it would have dropped 10–90%, and Yes below 10%. No ± values appear in this table.")
    A("")
    A("---")
    A("")

    # ---------------- Table S3
    A("## Table S3. Submitted (1,296-locus) versus co-processed (4,002-locus) values")
    A("")
    A("| Population | *n* | H_O submitted | H_O new | Δ H_O | Δ H_O (%) | H_E submitted | H_E new | "
      "Δ H_E | Δ H_E (%) | F_IS submitted | F_IS new | Polymorphic submitted | Polymorphic new | "
      "A_r(16) submitted | A_r(16) new | Δ A_r (%) |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    on = oldnew.set_index("population")
    for p in o:
        r = on.loc[p]
        A(f"| {PRETTY[p]} | {int(r.n)} | {f4(r.Ho_old)} | {f4(r.Ho_new)} | {s4(r.Ho_delta)} | "
          f"{r.Ho_pct_change:+.1f}% | {f4(r.He_old)} | {f4(r.He_new)} | {s4(r.He_delta)} | "
          f"{r.He_pct_change:+.1f}% | {s4(r.Fis_old)} | {s4(r.Fis_new)} | "
          f"{ithou(r.n_poly_old)} | {ithou(r.n_poly_new)} | {f4(r.Ar16_old)} | {f4(r.Ar16_new)} | "
          f"{'—' if pd.isna(r.Ar16_pct_change) else f'{r.Ar16_pct_change:+.1f}%'} |")
    A("")
    A("**Table S3.** Every population-level diversity value of the submitted manuscript beside "
      "its co-processed replacement. The submitted analysis used 1,296 loci recovered by "
      "intersecting CloneIDs between two separately processed DArT reports; the revised "
      f"analysis uses {ratios['loci_new']:,} loci from a single co-processed report "
      f"({ratios['locus_gain_factor']:.2f}×). "
      f"**Absolute H_E falls by {abs(ratios['mean_He_pct_change_mainland']):.1f}% on average "
      f"for the mainland populations and by {abs(ratios['He_pct_change_KI']):.1f}% for KI**, "
      "because the CloneID rule preferentially discarded low-diversity, near-fixed loci "
      "(mean H_E of the discarded loci 0.0581 against 0.1458 for those retained; Table S2). "
      "The relative conclusion is unchanged: KI H_E moves from "
      f"{ratios['KI_pct_of_mainland_mean_He_old']:.1f}% to "
      f"{ratios['KI_pct_of_mainland_mean_He_new']:.1f}% of the mainland mean and from "
      f"{ratios['KI_pct_of_OVENS_He_old']:.1f}% to {ratios['KI_pct_of_OVENS_He_new']:.1f}% of "
      "the Ovens value, and KI remains the lowest of the eleven populations. "
      "F_IS is not directly comparable between the two columns: the submitted value is the mean "
      "over loci of 1 − H_O/H_E, whereas the revised value is the Weir & Cockerham (1984) "
      "ratio-of-sums estimator. Dashes mark quantities not reported in the submitted "
      "manuscript. No ± values appear in this table.")
    A("")
    A("---")
    A("")

    # ---------------- Table S4
    A("## Table S4. MAF / MAC filtering sensitivity")
    A("")
    A("| Scenario | Loci | KI H_E | KI F_IS | Mainland mean H_E | Ovens H_E | "
      "KI as % of mainland mean | KI as % of Ovens | KI lowest of the eleven? |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---|")
    for scen, v in mafsens.items():
        A(f"| {scen} | {v['n_loci']:,} | {f4(v['by_pop']['KI']['He'])} | "
          f"{s4(v['by_pop']['KI']['Fis'])} | {f4(v['mainland_mean_He'])} | "
          f"{f4(v['by_pop']['OVENS']['He'])} | {pc1(v['KI_pct_of_mainland'])} | "
          f"{pc1(v['KI_pct_of_OVENS'])} | {'Yes' if v['KI_is_lowest_He'] else 'No'} |")
    A("")
    A("**Table S4.** Sensitivity of the headline comparison to minor-allele filtering, "
      "recomputed on the co-processed locus set. Thresholds were applied symmetrically across "
      "the joint 222-sample panel and never within KI alone. Absolute H_E rises steeply as rare "
      "variants are removed — from 0.0954 to 0.1899 for KI between no filter and MAF ≥ 0.05 — "
      "which is why absolute heterozygosities from differently filtered panels cannot be "
      "compared. The KI-to-mainland ratio is stable to within about seven percentage points "
      "across the whole range, and KI is the lowest of the eleven populations under every "
      "threshold. No ± values appear in this table.")
    A("")
    A("---")
    A("")

    # ---------------- Table S5
    A("## Table S5. Ahrens et al. (2025) Victorian benchmark")
    A("")
    A("| Ahrens population | *n* | H_O | H_E | N_e |")
    A("|---|---:|---:|---:|---:|")
    for name, n, ho, he, ne in [("Yarra", 220, 0.26, 0.28, 136.4), ("Tarago/Bunyip", 85, 0.26, 0.27, 53.8),
                                ("Maribyrnong", 71, 0.25, 0.26, 64.5), ("Olinda", 32, 0.23, 0.24, 8.6),
                                ("Lower Werribee", 28, 0.22, 0.22, 18.8), ("Upper Tarago", 25, 0.25, 0.26, 14.3),
                                ("Monbulk", 23, 0.22, 0.22, 5.4), ("Lang Lang", 11, 0.24, 0.24, 34.2),
                                ("Cardinia", 7, 0.21, 0.20, 18.3)]:
        A(f"| {name} | {n} | {ho:.2f} | {he:.2f} | {ne:.1f} |")
    A("")
    A("**Table S5.** Contemporary Victorian platypus diversity reported by Ahrens et al. (2025) "
      "for nine populations across five Melbourne-region catchments, used here only as an "
      "external, relative benchmark. Values are reproduced from that paper and are on the "
      "authors' own SNP scale (a Stacks-called, MAF- and linkage-filtered set); they are **not** "
      "comparable in absolute magnitude with the DArTseq-scale values in Table 1. The Ahrens "
      "dataset could not be co-processed with ours because its public deposit does not include "
      "raw reads. No ± values are reported by the source for these quantities.")
    A("")
    A("---")
    A("")

    # ---------------- Table S6
    A("## Table S6. Full 11 × 11 pairwise F_ST matrix")
    A("")
    pops = list(fstm.index)
    A("| | " + " | ".join(PRETTY[p] for p in pops) + " |")
    A("|---|" + "---:|" * len(pops))
    for a in pops:
        cells = []
        for b in pops:
            cells.append("—" if a == b else f3(fstm.loc[a, b]))
        A(f"| **{PRETTY[a]}** | " + " | ".join(cells) + " |")
    A("")
    A("**Pairwise F_ST with 95% bootstrap confidence intervals**")
    A("")
    A("| Pair | F_ST | 95% bootstrap CI |")
    A("|---|---:|---:|")
    for _, r in fstp.sort_values("Fst_WC").iterrows():
        A(f"| {PRETTY[r.pop_i]} vs {PRETTY[r.pop_j]} | {f3(r.Fst_WC)} | "
          f"[{f3(r.Fst_lo95)}, {f3(r.Fst_hi95)}] |")
    A("")
    A("**Table S6.** Complete matrix of pairwise Weir & Cockerham (1984) F_ST among the eleven "
      "populations on the 4,002-locus co-processed set, with all 55 pairwise values and their "
      "**1,000-replicate bootstrap-over-loci 95% confidence intervals** (no ± standard errors "
      "are reported). The mainland–mainland block is offered as a cross-panel check against "
      "Mijangos et al. (2022): unlike absolute heterozygosity, F_ST is a ratio of variance "
      "components and is far less sensitive to the ascertainment difference between the "
      "co-processed panel and the Mijangos-only panel. "
      f"Mean mainland–mainland F_ST is {f3(mm.Fst_WC.mean())} "
      f"(range {f3(mm.Fst_WC.min())}–{f3(mm.Fst_WC.max())}); mean KI–mainland F_ST is "
      f"{f3(kf.Fst_WC.mean())}.")
    A("")
    A("---")
    A("")

    # ---------------- Table S7
    A("## Table S7. KI pairwise kinship and the mainland nulls")
    A("")
    L_ = kin["KI"]["labels"]; Mx = kin["KI"]["matrix"]
    A("| | " + " | ".join(L_) + " |")
    A("|---|" + "---:|" * len(L_))
    for i, a in enumerate(L_):
        cells = ["—" if i == j else s3(Mx[i][j]) for j in range(len(L_))]
        A(f"| **{a}** | " + " | ".join(cells) + " |")
    A("")
    A("**Per-pair values with bootstrap confidence intervals and the null distributions**")
    A("")
    m8 = kin["matched_n8_null"]["kinship"]; fn = kin["empirical_null_n_ge_19"]
    A("| Distribution | *n* pairs | Mean | SD | 5th | 50th | 95th | 99th | Maximum |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    kv = np.array(kin["KI"]["values"], float)
    A(f"| **The 28 KI pairs** | 28 | {s3(kv.mean())} | {f3(kv.std(ddof=1))} | "
      f"{s3(np.percentile(kv,5))} | {s3(np.percentile(kv,50))} | {s3(np.percentile(kv,95))} | "
      f"{s3(np.percentile(kv,99))} | {s3(kv.max())} |")
    A(f"| Matched n = 8 mainland null (200 draws × 6 populations) | {m8['n']:,} | "
      f"{s3(m8['mean'])} | {f3(m8['sd'])} | {s3(m8['p05'])} | {s3(m8['p50'])} | "
      f"{s3(m8['p95'])} | {s3(m8['p99'])} | {s3(m8['max'])} |")
    A(f"| Full-sample mainland within-population null (n ≥ 19) | {fn['n_pairs']:,} | "
      f"{s3(fn['mean'])} | {f3(fn['sd'])} | {s3(fn['p05'])} | {s3(fn['p50'])} | "
      f"{s3(fn['p95'])} | {s3(fn['p99'])} | {s3(fn['max'])} |")
    A("")
    A("**Within-population kinship, per population (KI shown for comparison)**")
    A("")
    A("| Population | *n* | Pairs | Mean | SD | 5th | 50th | 95th |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|")
    for p, v in kin["within_pop_summary"].items():
        nm = "**Kangaroo Island**" if p == "KI" else PRETTY[p]
        A(f"| {nm} | {v['n']} | {v['n_pairs']:,} | {s3(v['mean'])} | {f3(v['sd'])} | "
          f"{s3(v['p05'])} | {s3(v['p50'])} | {s3(v['p95'])} |")
    A("")
    A("**Table S7.** Pairwise kinship among the eight Kangaroo Island platypuses and the null "
      "distributions against which it must be read. Values are Yang et al. (2010) genomic "
      "relationships on the **relatedness** scale (unrelated ≈ 0; full siblings and "
      "parent–offspring ≈ 0.5); halve for the kinship coefficient. KI values use KI allele "
      "frequencies. **The percentile columns are percentiles of the stated distribution, not "
      "± values; no standard errors or per-pair confidence intervals are reported in this "
      "table**. Per-pair bootstrap intervals were not carried into the analysis outputs, and a "
      "Mendelian simulation of full-sib, half-sib and parent–offspring pairs was not run; the "
      "matched n = 8 null below is therefore the basis on which any individual pair should be "
      "judged. Within-sample allele frequencies at n = 8 shift the "
      "whole distribution downward by construction, which is why the KI values must be compared "
      "with the matched n = 8 null rather than with a theoretical threshold. "
      f"{kin['KI']['n_above_matched_p95']} of the 28 KI pairs exceed the matched-n 95th "
      f"percentile and {kin['KI']['n_above_matched_max']} exceed its maximum; the empirical p "
      f"for the single highest KI pair is {kin['KI']['empirical_p_max_pair']:.4f}. "
      f"{kin['KI']['n_above_0p125_relatedness']} of 28 pairs exceed 0.125, the nominal "
      f"second-order threshold, and {kin['KI']['n_above_0p25_relatedness']} exceed 0.25. "
      "No KI pair reaches a relatedness consistent with first-order kinship.")
    A("")
    A("---")
    A("")

    # ---------------- Table S8
    A("## Table S8. Mainland per-individual metrics and the matched n = 8 nulls")
    A("")
    A("| Population | *n* | H_O mean | H_O range | sMLH mean | sMLH range | "
      "F_UNI (mainland frame) mean | F_UNI range |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|")
    for p in o:
        g = ind[ind.population == p]
        sm = g.Ho_ind / panel_mean_ho
        nm = "**Kangaroo Island**" if p == "KI" else PRETTY[p]
        A(f"| {nm} | {len(g)} | {f4(g.Ho_ind.mean())} | {f4(g.Ho_ind.min())}–{f4(g.Ho_ind.max())} | "
          f"{f3(sm.mean())} | {f3(sm.min())}–{f3(sm.max())} | {s3(g.F_uni_mainland.mean())} | "
          f"{s3(g.F_uni_mainland.min())} to {s3(g.F_uni_mainland.max())} |")
    A(f"| *All mainland (214)* | 214 | {f4(mi.Ho_ind.mean())} | "
      f"{f4(mi.Ho_ind.min())}–{f4(mi.Ho_ind.max())} | {f3((mi.Ho_ind/panel_mean_ho).mean())} | "
      f"{f3((mi.Ho_ind/panel_mean_ho).min())}–{f3((mi.Ho_ind/panel_mean_ho).max())} | "
      f"{s3(mi.F_uni_mainland.mean())} | {s3(mi.F_uni_mainland.min())} to "
      f"{s3(mi.F_uni_mainland.max())} |")
    A("")
    A("**Matched n = 8 null distributions for individual inbreeding (own-population frame)**")
    A("")
    fu = kin["matched_n8_null"]["F_uni"]; fh = kin["matched_n8_null"]["F_hom"]
    fg = kin["matched_n8_null"]["F_grm"]
    A("| Estimator | *n* values | Mean | SD | 5th percentile | 95th percentile | Maximum |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    A(f"| F_UNI (Yang et al. 2010, F-hat-3) | {fu['n']:,} | {s3(fu['mean'])} | {f3(fu['sd'])} | "
      f"{s3(fu['p05'])} | {s3(fu['p95'])} | {s3(fu['max'])} |")
    A(f"| F_HOM (excess homozygosity) | {fu['n']:,} | {s3(fh['mean'])} | — | {s3(fh['p05'])} | "
      f"{s3(fh['p95'])} | — |")
    A(f"| F_GRM (F-hat-1) | {fu['n']:,} | {s3(fg['mean'])} | — | {s3(fg['p05'])} | "
      f"{s3(fg['p95'])} | — |")
    A("")
    A("**Concordance of the three individual-F estimators for KI1–KI8**")
    A("")
    A("| Individual | F_HOM (KI frame) | F_UNI (KI frame) | F_GRM (KI frame) | "
      "F_HOM (mainland frame) | F_UNI (mainland frame) | F_GRM (mainland frame) |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    for _, r in kid.iterrows():
        A(f"| {r['sample']} | {s3(r.F_hom_ownpop)} | {s3(r.F_uni_ownpop)} | {s3(r.F_grm_ownpop)} | "
          f"{s3(r.F_hom_mainland)} | {s3(r.F_uni_mainland)} | {s3(r.F_grm_mainland)} |")
    A("")
    A("**Table S8.** Per-individual heterozygosity and inbreeding summarised by population, the "
      "matched n = 8 null distributions for individual F, and the concordance of the three "
      "individual-F estimators for the eight KI animals. sMLH is each individual's H_O divided "
      f"by the mean H_O of all 222 individuals ({f4(panel_mean_ho)}). **Ranges are minima and "
      "maxima and the percentile columns are percentiles of the null distribution; no ± values "
      "appear in this table.** The nulls were built from 200 draws of n = 8 from each of the six "
      "populations with n ≥ 19. F_GRM is the most rare-allele-sensitive of the three estimators "
      "and inflates strongly when evaluated in a foreign reference frame, so mainland-frame "
      "F_GRM values should be read as a rank ordering rather than as absolute inbreeding "
      "coefficients; F_HOM and F_UNI are the estimators quoted in the main text. "
      f"Per-individual H_O ranges from {f4(ind[ind.population=='KI'].Ho_ind.min())} to "
      f"{f4(ind[ind.population=='KI'].Ho_ind.max())} in KI and from {f4(mi.Ho_ind.min())} to "
      f"{f4(mi.Ho_ind.max())} on the mainland: every KI animal sits below the mainland mean and "
      f"the KI mean falls at the {jsum['individuals']['KI_mean_percentile_in_mainland']:.1f}th "
      f"percentile of the mainland distribution, but the ranges overlap — "
      f"{jsum['individuals']['n_mainland_below_KI_max']} of the 214 mainland animals have a lower "
      f"H_O than the highest KI animal and "
      f"{jsum['individuals']['n_KI_below_all_mainland_Ho']} KI animals fall below the mainland "
      "minimum. The KI deficit is a shift in the whole distribution, not a set of individually "
      "exceptional animals.")
    A("")
    A("---")
    A("")

    # ---------------- Table S9
    sx = cascade["sexchrom"]; cm = sx["confusion_vs_direct"]
    A("## Table S9. Scaffold→chromosome crosswalk and sex-linkage validation")
    A("")
    A("| Quantity | Value |")
    A("|---|---:|")
    A(f"| Method | {sx['method'].replace('_',' ')} |")
    A(f"| AlleleID anchors shared with the Mijangos-only report | {sx['n_anchor_loci']:,} |")
    A(f"| PTTO01 scaffolds in the co-processed report | {sx['n_scaffolds_total']:,} |")
    A(f"| Scaffolds anchored to a RefSeq chromosome | {sx['n_scaffolds_anchored']:,} |")
    A(f"| Anchored scaffolds ≥ 95% pure for one chromosome | {sx['n_scaffolds_pure_ge95pct']:,} |")
    A(f"| Locus-level agreement with direct annotation | {100*sx['validation_agreement_vs_direct']:.2f}% |")
    A(f"| — direct sex-linked, crosswalk sex-linked | {cm['direct_sex_cw_sex']:,} |")
    A(f"| — direct sex-linked, crosswalk autosomal | {cm['direct_sex_cw_auto']:,} |")
    A(f"| — direct autosomal, crosswalk sex-linked | {cm['direct_auto_cw_sex']:,} |")
    A(f"| — direct autosomal, crosswalk autosomal | {cm['direct_auto_cw_auto']:,} |")
    A(f"| Sex-linked loci removed | {sx['n_sex_loci']:,} (11.7% of the input panel) |")
    A(f"| Loci on scaffolds with no crosswalk anchor | {sx['unanchored_loci']:,} |")
    A(f"| — final loci if these are dropped | {sx['sensitivity_drop_unanchored']['final_loci']:,} |")
    A(f"| — largest resulting shift in any population H_E | {sx['sensitivity_drop_unanchored']['max_He_shift']:.5f} |")
    A(f"| Loci flagged by the data-driven male/female heterozygosity test | {sx['data_driven_candidates']:,} |")
    A(f"| — of these, on autosomal scaffolds | {sx['data_driven_candidates_on_autosomal_scaffolds']:,} |")
    A(f"| — of these, surviving into the final set | {sx['data_driven_candidates_in_final_set']:,} |")
    A("")
    A("**Table S9.** Construction and validation of the scaffold→chromosome crosswalk required "
      "for sex-chromosome exclusion. The co-processed DArT report maps tags to PTTO01 scaffolds, "
      "whereas the Mijangos-only report used for the submitted analysis mapped them to RefSeq "
      "chromosomes, so a crosswalk was built from the AlleleID anchors the two reports share and "
      "each scaffold assigned to its majority chromosome. It was validated three ways: scaffold "
      "purity, locus-level agreement against the older report's direct annotation, and a "
      "data-driven male/female heterozygosity test. Residual sex-linkage candidates are reported "
      "rather than filtered on: only eight survive into the final 4,002-locus set. **No ± values "
      "appear in this table**; the single sensitivity figure quoted is the largest absolute shift "
      "in any population's H_E when the 347 unanchored loci are dropped.")
    A("")
    A("---")
    A("")
    A("*Symbols: H_O observed heterozygosity; H_E expected heterozygosity; F_IS inbreeding "
      "coefficient within populations; F_ST differentiation among populations; N_e effective "
      "population size; A_r rarefied allelic richness; sMLH standardised multilocus "
      "heterozygosity. Subscripts are to be set as true subscripts in the typeset manuscript.*")
    A("")

    (WORK / "TABLES.md").write_text("\n".join(L), encoding="utf-8")
    print(f"  wrote TABLES.md ({(WORK/'TABLES.md').stat().st_size:,} B)")


# ======================================================================= main
if __name__ == "__main__":
    ov = ovens_full_spectrum()
    figure1()
    figure2()
    figure3()
    figure4(ov)
    figureS1()
    figureS2()
    figureS3(ov)
    tables()
    print("done.")
