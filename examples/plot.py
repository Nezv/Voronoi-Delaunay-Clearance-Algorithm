"""Draws each example: the nets, the curve shaded by how much of its requirement each gap uses, and every violation
as its gap (red). Needs matplotlib. Writes docs/examples/<example>.png."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from shapely.geometry import MultiPolygon, box

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import facing_edges
import levels
import openings
import wavy_edge
from clearance_curve import ClearanceCurve

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "examples")
COLOURS = ["#f1d3c7", "#cfdcf3", "#eee0b8", "#d5ead9"]
EDGES = ["#b8432d", "#2c5bc7", "#8a6516", "#2e7d45"]


def draw(name, nets, requirement, bare_requirement=None, openings_=None, step=0.02, tau=0.02, view=None):
    curve = ClearanceCurve(nets, step=step, openings=openings_)
    found = curve.check(requirement, bare_requirement, tau=tau)
    fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
    if openings_ is not None:
        ax.fill(*openings_.exterior.xy, color="#e9e9e9", zorder=0, label="opening")
    for k, (net, geom) in enumerate(nets.items()):
        for poly in (geom.geoms if isinstance(geom, MultiPolygon) else [geom]):
            ax.fill(*poly.exterior.xy, color=COLOURS[k % 4], ec=EDGES[k % 4], lw=1, zorder=1)
        c = (geom.intersection(box(*view)) if view else geom).representative_point()
        ax.text(c.x, c.y, net, ha="center", va="center", fontsize=9, color="#444", zorder=4)
    seg, g = curve.edges()
    use = curve.ratio(requirement, bare_requirement)
    seg, use = seg[g < 1.0], use[g < 1.0]                 # the curve where nets are under 1 mm apart
    lc = LineCollection(seg, cmap="Greys", norm=plt.Normalize(0.3, 1.2), linewidths=1.6, zorder=2)
    lc.set_array(np.minimum(use, 1.2))
    ax.add_collection(lc)
    for v in found:
        ax.plot([v.foot_a[0], v.foot_b[0]], [v.foot_a[1], v.foot_b[1]], color="#d62728", lw=2, zorder=3)
    ax.set_aspect("equal")
    if view:
        ax.set_xlim(view[0], view[2])
        ax.set_ylim(view[1], view[3])
    ax.set_title(f"{name}: {len(found)} violation{'s' if len(found) != 1 else ''}", fontsize=10)
    ax.set_xlabel("mm")
    fig.colorbar(lc, ax=ax, shrink=0.7, label="requirement / gap")
    fig.tight_layout()
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, f"{name}.png"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    draw("facing_edges", facing_edges.nets, 0.2, view=(-0.2, -0.2, 5.2, 2.7))
    draw("wavy_edge", wavy_edge.nets, 0.2, step=0.005, tau=0.03, view=(-1, -0.3, 1, 0.5))
    draw("levels", levels.nets, levels.requirement, view=(-1.5, -0.8, 1.5, 0.8))
    draw("openings", openings.nets, 0.13, 0.25, openings.openings, view=(-0.2, -0.4, 6.4, 1.4))
    print("written to", os.path.normpath(OUT))
