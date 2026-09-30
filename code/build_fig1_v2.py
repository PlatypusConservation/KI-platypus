#!/usr/bin/env python3
"""
build_fig1_v2.py
================
Rebuilds Figure 1 (two-panel study-area map) of the Kangaroo Island platypus
conservation-genomics manuscript in the house style used by
build_revision_figures_v2.py for Figures 2-4.

Panel (a)  south-eastern Australia: the Kangaroo Island recipient population
           (plotted at the mean of the panel-(b) capture sites, i.e. on the
           Rocky River at the island's western end),
           the two founder sources (Healesville / upper Yarra, Vic., 16
           founders; Wynyard, Tas., 3 founders) and the mainland
           reference-population centroids grouped by region.
Panel (b)  the Rocky River, Flinders Chase National Park, with the four sites
           at which the eight genotyped individuals were captured; marker area
           is proportional to the number of individuals caught at the site.

House style (differences from the original build_fig1_final.py):
  * no titles, no sub-headings, no year text anywhere in the figure;
  * no text of any kind inside the plotting areas (the per-site "x4/x2/x1/x1"
    annotations are replaced by proportional marker areas) - only the scale
    bars carry their own labels;
  * a single horizontal legend for the whole figure, below both axes;
  * bare "(a)"/"(b)" panel identifiers in the outside top-left corner;
  * Okabe-Ito colourblind-safe palette and marker shapes shared with
    Figures 2-4: Border Rivers = blue circle, Snowy = green square,
    Upper Murray = pink triangle, Kangaroo Island = vermillion diamond.

Coastline: GIS DATA/swma_geom_nat_polygon.shp (the original's
australia_bndGDA94.shp is no longer available).  Polygons are filled with the
land colour and stroked with the same colour, so only the national silhouette
shows and no internal catchment boundaries are visible.

Python 3 + numpy + pyshp + matplotlib only.
"""

import csv
import collections
import os
from pathlib import Path

import numpy as np
import shapefile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.lines import Line2D

# --------------------------------------------------------------------- paths
UP   = Path("/mnt/user-data/uploads")
GIS  = UP / "GIS DATA"
KI   = UP / "Kangaroo Island"
GEN  = KI / "Genetics"
OUT  = Path("/mnt/user-data/outputs")
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

COAST = GIS / "swma_geom_nat_polygon"      # .shp + .shx, no .dbf needed

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

# Okabe-Ito, identical to Figures 2-4
C = {"BORDER": "#0072B2", "SNOWY": "#009E73", "U_MURRAY": "#CC79A7", "KI": "#D55E00"}
M = {"BORDER": "o",       "SNOWY": "s",       "U_MURRAY": "^",       "KI": "D"}
REGION_LABEL = {"BORDER":   "Border Rivers (NSW)",
                "SNOWY":    "Snowy (NSW)",
                "U_MURRAY": "Upper Murray (Vic.)"}
VIC_SRC = "#F0E442"     # Okabe-Ito yellow  - Healesville / upper Yarra founders
TAS_SRC = "#56B4E9"     # Okabe-Ito sky blue - Wynyard founders
LAND    = "#e9e4d8"
SEA     = "#ffffff"
RIVER   = "#0072B2"
GREY    = "#4d4d4d"

# Founder-source localities (town centroids, approximate - as in the original)
HEALESVILLE = (145.5172, -37.6540)
WYNYARD     = (145.7262, -40.9897)

# ------------------------------------------------- reference-population centroids
ind2pop = {r[0].strip(): r[1].strip()
           for r in csv.reader(open(GEN / "mijangos_data/new_pop_assignments.csv"))
           if len(r) >= 2}
pop2reg = {r[0].strip(): r[1].strip()
           for r in csv.reader(open(GEN / "mijangos_data/new_pop_assignments_regions.csv"))
           if len(r) >= 2}
coords = collections.defaultdict(list)
for row in csv.DictReader(open(GEN / "mijangos_data/ID_Pop_Platypus.csv")):
    p = ind2pop.get(row["pop"].strip())
    try:
        lat, lon = float(row["lat"]), float(row["lon"])
    except (TypeError, ValueError):
        continue
    if p:
        coords[p].append((lon, lat))
cent = {p: (float(np.mean([c[0] for c in v])), float(np.mean([c[1] for c in v])))
        for p, v in coords.items()}

# ------------------------------------------------ Rocky River capture sites
rr = shapefile.Reader(str(KI / "GIS/RockyRiver.shp"))

sv = shapefile.Reader(str(KI / "GIS/Platypus_Survey_Sites.shp"))
fld = [f[0] for f in sv.fields[1:]]
cap = [dict(zip(fld, list(rec))) for rec in sv.iterRecords()]
cap = [c for c in cap if int(c.get("Platypus21") or 0) > 0]
cap.sort(key=lambda c: -int(c["Platypus21"]))
cap_lon = [c["Long"] for c in cap]
cap_lat = [c["Lat"] for c in cap]
cap_n   = [int(c["Platypus21"]) for c in cap]
assert sum(cap_n) == 8 and len(cap_n) == 4, (cap_n,)

# The Kangaroo Island recipient population is marked in panel (a) at the mean
# of those capture sites - the Rocky River, Flinders Chase, at the western end
# of the island - rather than at a hard-coded island centroid.
KI_POP = (float(np.mean(cap_lon)), float(np.mean(cap_lat)))


# ------------------------------------------------------------------- geometry
def poly_parts(shape):
    pts = shape.points
    pp = list(shape.parts) + [len(pts)]
    return [pts[pp[i]:pp[i + 1]] for i in range(len(pp) - 1)]


def land_rings(bbox):
    """All coastline polygon rings whose shape bbox intersects `bbox`."""
    x0, y0, x1, y1 = bbox
    sf = shapefile.Reader(shp=open(str(COAST) + ".shp", "rb"),
                          shx=open(str(COAST) + ".shx", "rb"))
    rings = []
    for sh in sf.iterShapes():
        b = sh.bbox
        if b[2] < x0 or b[0] > x1 or b[3] < y0 or b[1] > y1:
            continue
        rings.extend(poly_parts(sh))
    return rings


def simplify(ring, tol):
    """Drop consecutive vertices closer than `tol` degrees.  At the scale of
    panel (a) this is invisible but keeps the vector PDF small."""
    if tol <= 0 or len(ring) < 4:
        return ring
    out = [ring[0]]
    for pt in ring[1:-1]:
        if abs(pt[0] - out[-1][0]) + abs(pt[1] - out[-1][1]) >= tol:
            out.append(pt)
    out.append(ring[-1])
    return out if len(out) >= 3 else ring


def point_in_rings(pt, rings):
    """True if `pt` falls inside any of `rings` (even-odd ray casting)."""
    x, y = pt
    for ring in rings:
        inside = False
        n = len(ring)
        for i in range(n):
            x0, y0 = ring[i]
            x1, y1 = ring[(i + 1) % n]
            if (y0 > y) != (y1 > y):
                xi = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
                if x < xi:
                    inside = not inside
        if inside:
            return True
    return False


def draw_land(ax, bbox, lw=0.6, tol=0.0):
    """Fill the land silhouette; stroke in the land colour so that the seams
    between adjacent catchment polygons close and no internal edges show."""
    rings = [simplify(r, tol) for r in land_rings(bbox)]
    ax.add_collection(PolyCollection(rings, facecolor=LAND, edgecolor=LAND,
                                     linewidths=lw, zorder=1))
    return rings


def scale_bar(ax, km, lat_ref, frac_x=0.055, frac_y=0.075, lw=2.2):
    """Horizontal scale bar with its own label (the only text allowed inside
    a panel under the house style)."""
    deg = km / (111.32 * np.cos(np.deg2rad(lat_ref)))
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    x = x0 + (x1 - x0) * frac_x
    y = y0 + (y1 - y0) * frac_y
    ax.plot([x, x + deg], [y, y], color="black", lw=lw, zorder=9,
            solid_capstyle="butt")
    ax.text(x + deg / 2, y + (y1 - y0) * 0.022, f"{km:g} km",
            fontsize=7.6, ha="center", va="bottom", zorder=9)


def panel_labels(fig, axes, letters):
    """Bare parenthesised identifiers in the conventional outside top-left
    corner, drawn in figure coordinates (so nothing is clipped) and on a
    common left edge so that the two panels agree."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bbs = [ax.get_tightbbox(r).transformed(fig.transFigure.inverted())
           for ax in axes]
    x = min(bb.x0 for bb in bbs) - 0.012
    for bb, letter in zip(bbs, letters):
        fig.text(x, bb.y1 + 0.005, f"({letter})", fontsize=10.5,
                 fontweight="bold", ha="left", va="bottom")


# =============================================================== the figure
fig = plt.figure(figsize=(7.0, 8.85))
gs = fig.add_gridspec(2, 1, height_ratios=[1.34, 1.0], hspace=0.22,
                      left=0.11, right=0.985, top=0.945, bottom=0.155)

# ------------------------------------------------------------------ panel (a)
axA = fig.add_subplot(gs[0])
axA.set_facecolor(SEA)
XA0, XA1, YA0, YA1 = 133.5, 154.5, -44.5, -27.8
draw_land(axA, (XA0 - 1, YA0 - 1, XA1 + 1, YA1 + 1), lw=0.6, tol=0.006)

# the recipient-population marker must fall on Kangaroo Island, not in the sea
assert point_in_rings(KI_POP, land_rings((136.3, -36.3, 138.3, -35.5))), KI_POP

for p, (lon, lat) in cent.items():
    reg = pop2reg.get(p)
    if reg not in C:
        continue
    axA.scatter([lon], [lat], s=16, marker=M[reg], facecolor=C[reg],
                edgecolor="white", linewidths=0.4, zorder=5)

axA.scatter([HEALESVILLE[0]], [HEALESVILLE[1]], s=60, marker="*",
            facecolor=VIC_SRC, edgecolor="black", linewidths=0.5, zorder=6)
axA.scatter([WYNYARD[0]], [WYNYARD[1]], s=60, marker="*",
            facecolor=TAS_SRC, edgecolor="black", linewidths=0.5, zorder=6)
axA.scatter([KI_POP[0]], [KI_POP[1]], s=22, marker="D", facecolor=C["KI"],
            edgecolor="black", linewidths=0.45, zorder=7)

axA.set_xlim(XA0, XA1)
axA.set_ylim(YA0, YA1)
axA.set_aspect(1 / np.cos(np.deg2rad(36)))
axA.set_xlabel("Longitude (°E)")
axA.set_ylabel("Latitude (°S)")
axA.set_yticks([-44, -42, -40, -38, -36, -34, -32, -30, -28])
axA.set_yticklabels([f"{abs(t):g}" for t in axA.get_yticks()])
axA.set_xticks([134, 138, 142, 146, 150, 154])
scale_bar(axA, 200, 36.0)

# ------------------------------------------------------------------ panel (b)
axB = fig.add_subplot(gs[1])

segs = [part for sh in rr.iterShapes() for part in poly_parts(sh)]
allx = cap_lon + [p[0] for s in segs for p in s]
ally = cap_lat + [p[1] for s in segs for p in s]
padx = (max(allx) - min(allx)) * 0.07
pady = (max(ally) - min(ally)) * 0.12
extra_s = (max(ally) - min(ally)) * 0.25      # show the coastline at the SW
XB0, XB1 = min(allx) - padx, max(allx) + padx
YB0, YB1 = min(ally) - pady - extra_s, max(ally) + pady

axB.set_facecolor("#dfe9f0")                  # sea
draw_land(axB, (XB0 - 0.05, YB0 - 0.05, XB1 + 0.05, YB1 + 0.05), lw=0.3)
axB.add_collection(LineCollection(segs, colors=RIVER, linewidths=1.4, zorder=3))

SIZE_PER_ANIMAL = 40.0                        # marker area proportional to n
axB.scatter(cap_lon, cap_lat, s=[SIZE_PER_ANIMAL * n for n in cap_n],
            marker="D", facecolor=C["KI"], edgecolor="black", linewidths=0.7,
            zorder=6)

axB.set_xlim(XB0, XB1)
axB.set_ylim(YB0, YB1)
axB.set_aspect(1 / np.cos(np.deg2rad(35.95)))
axB.set_xlabel("Longitude (°E)")
axB.set_ylabel("Latitude (°S)")
axB.set_xticks([136.66, 136.70, 136.74, 136.78, 136.82])
axB.set_yticks([-35.98, -35.96, -35.94, -35.92, -35.90, -35.88])
axB.set_yticklabels([f"{abs(t):.2f}" for t in axB.get_yticks()])
scale_bar(axB, 5, 35.95)

# --------------------------------------------- one legend, horizontal, at foot
handles = [
    Line2D([0], [0], marker="D", color="none", markerfacecolor=C["KI"],
           markeredgecolor="black", markersize=7.0,
           label="Kangaroo Island recipient population"),
    Line2D([0], [0], marker="*", color="none", markerfacecolor=VIC_SRC,
           markeredgecolor="black", markersize=12.0,
           label="Healesville / upper Yarra (Vic.), 16 founders"),
    Line2D([0], [0], marker="*", color="none", markerfacecolor=TAS_SRC,
           markeredgecolor="black", markersize=12.0,
           label="Wynyard (Tas.), 3 founders"),
]
for r in ("BORDER", "SNOWY", "U_MURRAY"):
    handles.append(Line2D([0], [0], marker=M[r], color="none",
                          markerfacecolor=C[r], markeredgecolor="white",
                          markersize=6.0,
                          label=REGION_LABEL[r] + " reference populations"))
handles.append(Line2D([0], [0], color=RIVER, lw=1.6, label="Rocky River"))
for n in (1, 2, 4):
    handles.append(Line2D([0], [0], marker="D", color="none",
                          markerfacecolor=C["KI"], markeredgecolor="black",
                          markersize=np.sqrt(SIZE_PER_ANIMAL * n) * 0.80,
                          label=("Capture site, %d individual%s"
                                 % (n, "" if n == 1 else "s"))))

panel_labels(fig, (axA, axB), ("a", "b"))

leg = fig.legend(handles=handles, loc="lower center", ncol=3,
                 bbox_to_anchor=(0.5, 0.006), handletextpad=0.6,
                 columnspacing=1.4, labelspacing=0.65, borderaxespad=0.0)

for stem in ("Figure_1_study_area",):
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches=None, facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches=None, facecolor="white")
plt.close(fig)

print("capture sites (Platypus21 > 0):")
for c, n in zip(cap, cap_n):
    print(f"  {c['Site']:>5}  {c['SiteName']:<9}  lat {c['Lat']:.6f}  "
          f"lon {c['Long']:.6f}  n = {n}")
print(f"  total {sum(cap_n)} individuals across {len(cap_n)} sites")
print("wrote Figure_1_study_area.png / .pdf")
