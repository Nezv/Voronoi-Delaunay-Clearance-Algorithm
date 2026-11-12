"""The clearance curve: where two conductors are equally near, and the gap across it.

Every net's outline is sampled every `step` mm and the samples are triangulated (Delaunay). A triangle edge whose
two ends lie on different nets crosses the curve: its dual edge in the Voronoi diagram of the samples is a piece of
the curve, and its two ends are the foot points of a clear straight gap between those nets. The triangulation is
built once per geometry; checks against any requirement re-use it.

Units are whatever the geometry is in; the defaults (step, tau, eps) assume millimetres.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Union

import numpy as np
import shapely
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import Delaunay, cKDTree

Requirement = Union[float, Callable[[str, str], float]]

EPS = 0.002              # a gap within this of its requirement passes
LEVEL_RADIUS = 0.2       # level-curve chords of one pair closer than this are one violation
OPENING_GROW = 0.005     # openings are grown by this so an edge on their border counts as bare


@dataclass(frozen=True)
class Violation:
    """One connected stretch of curve under its requirement, given at its worst point."""
    net_a: str
    net_b: str
    gap: float                 # the gap at the worst point
    required: float            # its requirement
    foot_a: tuple              # the gap's ends, on net_a and net_b
    foot_b: tuple
    method: str                # "curve": on the curve of all nets; "level": found by the level search (section 5)
    length: float              # length of curve in the violation (0 for "level")
    box: tuple                 # (x0, y0, x1, y1) around it
    bare: bool = False         # a foot point lies in an opening

    @property
    def ratio(self):
        """U = c / g: above 1 is a violation."""
        return self.required / self.gap if self.gap > 0 else float("inf")


@dataclass(frozen=True)
class Gap:
    """The narrowest gap of one stretch of curve between two nets."""
    net_a: str
    net_b: str
    gap: float
    foot_a: tuple
    foot_b: tuple
    bare: bool = False


class ClearanceCurve:
    """The clearance curve of a set of nets on one layer.

    nets      {name: shapely geometry}, the copper of each net (polygons; lines and points have no area but work)
    step      sampling distance along every outline
    openings  optional geometry: where the copper is bare (for example solder mask openings). Gaps with a foot
              point inside it are checked against bare_requirement.
    max_gap   gaps longer than this are dropped when the curve is built: the largest requirement ever checked
    """

    def __init__(self, nets: Mapping[str, object], step: float = 0.02, openings=None, max_gap: float = 10.0):
        self.step = float(step)
        self.max_gap = float(max_gap)
        self.nets = [n for n, g in nets.items() if g is not None and not g.is_empty]
        self.geoms = [nets[n] for n in self.nets]
        pts, lab = [], []
        for i, geom in enumerate(self.geoms):
            c = shapely.get_coordinates(shapely.segmentize(_outline(geom), self.step))
            pts.append(c)
            lab.append(np.full(len(c), i, np.int32))
        P = np.vstack(pts) if pts else np.zeros((0, 2))
        L = np.concatenate(lab) if lab else np.zeros(0, np.int32)
        _, first = np.unique(np.round(P * 1e5).astype(np.int64), axis=0, return_index=True)
        first.sort()
        self.P, self.L = P[first], L[first]
        if openings is not None and not openings.is_empty and len(self.P):
            self.bare = shapely.contains_xy(openings.buffer(OPENING_GROW), self.P[:, 0], self.P[:, 1])
        else:
            self.bare = np.zeros(len(self.P), bool)
        self._ridges()

    # ------------------------------------------------------------------------------------------ the curve
    def _ridges(self):
        P, L = self.P, self.L
        empty = np.zeros(0, int)
        self.T = self.nb = self.a = self.b = self.la = self.lb = empty
        self.g = np.zeros(0)
        self.V1 = self.V2 = np.zeros((0, 2))
        if len(P) < 3 or len(self.nets) < 2:
            return
        try:
            tri = Delaunay(P - P.mean(0))
        except Exception:                            # every sample on one line: no curve
            return
        S, NB = tri.simplices, tri.neighbors
        A, B, C = P[S[:, 0]], P[S[:, 1]], P[S[:, 2]]
        d = 2 * (A[:, 0] * (B[:, 1] - C[:, 1]) + B[:, 0] * (C[:, 1] - A[:, 1]) + C[:, 0] * (A[:, 1] - B[:, 1]))
        a2, b2, c2 = (A ** 2).sum(1), (B ** 2).sum(1), (C ** 2).sum(1)
        with np.errstate(divide="ignore", invalid="ignore"):       # circumcentres: the Voronoi vertices
            CC = np.stack([(a2 * (B[:, 1] - C[:, 1]) + b2 * (C[:, 1] - A[:, 1]) + c2 * (A[:, 1] - B[:, 1])) / d,
                           (a2 * (C[:, 0] - B[:, 0]) + b2 * (A[:, 0] - C[:, 0]) + c2 * (B[:, 0] - A[:, 0])) / d], 1)
        M = len(S)
        T = np.repeat(np.arange(M), 3)
        I = np.tile(np.arange(3), M)
        nb = NB.ravel()
        k = nb > T                                   # each inner edge once
        T, I, nb = T[k], I[k], nb[k]
        a, b = S[T, (I + 1) % 3], S[T, (I + 2) % 3]
        la, lb = L[a], L[b]
        g = np.linalg.norm(P[a] - P[b], axis=1)
        k = (la != lb) & (g < self.max_gap) & np.isfinite(CC[T]).all(1) & np.isfinite(CC[nb]).all(1)
        self.T, self.nb, self.a, self.b, self.la, self.lb, self.g = T[k], nb[k], a[k], b[k], la[k], lb[k], g[k]
        self.V1, self.V2 = CC[self.T], CC[self.nb]

    def edges(self):
        """The curve as an (n, 2, 2) array of line segments, and the gap across each: for drawing it."""
        return np.stack([self.V1, self.V2], 1), self.g.copy()

    def ratio(self, requirement: Requirement, bare_requirement: Requirement = None):
        """U = c / g along the curve (section 5), one value per segment of edges(): above 1 is under the
        requirement. Hidden pairs (the level search) are not on the curve and not in this."""
        C, Cb = self._matrices(requirement, bare_requirement)
        bare = self.bare[self.a] | self.bare[self.b]
        c = np.where(bare, Cb[self.la, self.lb], C[self.la, self.lb])
        return c / np.maximum(self.g, 1e-12)

    # ------------------------------------------------------------------------------------------ measuring
    def gaps(self, threshold: float):
        """Every gap under the threshold: one per stretch of curve between two nets, at its narrowest point.
        Narrowest first."""
        k = np.nonzero(self.g < threshold)[0]
        if not len(k):
            return []
        T, nb, la, lb, g = self.T[k], self.nb[k], self.la[k], self.lb[k], self.g[k]
        n, N = len(k), len(self.nets)
        pair = np.minimum(la, lb).astype(np.int64) * (N + 1) + np.maximum(la, lb)
        M = int(max(T.max(), nb.max())) + 1
        nodes, inv = np.unique(np.concatenate([pair * M + T, pair * M + nb]), return_inverse=True)
        inv = inv.ravel()
        gr = coo_matrix((np.ones(n), (inv[:n], inv[n:])), shape=(len(nodes), len(nodes)))
        _, comp = connected_components(gr, directed=False)
        rc = comp[inv[:n]]
        order = np.lexsort((g, rc))
        first = order[np.r_[True, rc[order][1:] != rc[order][:-1]]]
        out = []
        for w in first:
            i, j, fa, fb = self._ordered(int(la[w]), int(lb[w]), self.P[self.a[k[w]]], self.P[self.b[k[w]]])
            out.append(Gap(self.nets[i], self.nets[j], float(g[w]), fa, fb,
                           bool(self.bare[self.a[k[w]]] or self.bare[self.b[k[w]]])))
        return sorted(out, key=lambda x: x.gap)

    # ------------------------------------------------------------------------------------------ checking
    def check(self, requirement: Requirement, bare_requirement: Requirement = None, tau: float = 0.02,
              eps: float = EPS):
        """Every violation of the requirement, worst first.

        requirement       one value for every pair, or a function (net_a, net_b) -> value
        bare_requirement  the same, for gaps with a foot point in the openings (default: requirement)
        tau               the band: dips of one stretch of curve less than tau above the requirement apart are
                          one violation (section 6)
        eps               a gap within eps of its requirement passes
        """
        C, Cb = self._matrices(requirement, bare_requirement)
        if not len(self.g):
            return []
        cmax = float(max(C.max(), Cb.max()))
        k = self.g < cmax + tau
        T, nb, a, b, la, lb, g, V1, V2 = (x[k] for x in (self.T, self.nb, self.a, self.b, self.la, self.lb, self.g,
                                                         self.V1, self.V2))
        bare = self.bare[a] | self.bare[b]
        c = np.where(bare, Cb[la, lb], C[la, lb])
        viol = g < c - eps
        out = []
        for idx in _zones(T, nb, g < c + tau, viol):
            w = idx[np.argmax(c[idx] / np.maximum(g[idx], 1e-12))]
            fa, fb = self.P[a[w]], self.P[b[w]]
            seg = np.vstack([V1[idx], V2[idx], fa, fb])
            box = (seg[:, 0].min(), seg[:, 1].min(), seg[:, 0].max(), seg[:, 1].max())
            length = float(np.linalg.norm(V1[idx] - V2[idx], axis=1).sum())
            out.append(self._violation("curve", int(la[w]), int(lb[w]), g[w], c[w], fa, fb, box, length, bare[w]))
        seen = (la[viol], lb[viol], a[viol], b[viol])
        out += self._levels(C, Cb, seen, eps)
        return sorted(out, key=lambda v: -v.ratio)

    def _matrices(self, requirement, bare_requirement):
        C = self._matrix(requirement)
        return C, (C if bare_requirement is None else self._matrix(bare_requirement))

    def _matrix(self, req):
        """The requirement between every two nets, as a matrix (0 on the diagonal)."""
        N = len(self.nets)
        if callable(req):
            C = np.zeros((N, N))
            for i in range(N):
                for j in range(i + 1, N):
                    C[i, j] = C[j, i] = float(req(self.nets[i], self.nets[j]))
        else:
            C = np.full((N, N), float(req))
            np.fill_diagonal(C, 0.0)
        return C

    def _levels(self, C, Cb, seen, eps):
        """Section 5: pairs a net with a smaller requirement hides from the curve. Each sample of a net that needs
        more than the smallest requirement looks for other nets within its largest requirement; the nearest
        violating partner per sample and net is kept when the curve doesn't already show that pair nearby and the
        straight gap crosses no other copper."""
        P, L, X = self.P, self.L, self.bare
        N = len(self.nets)
        if N < 2:
            return []
        off = ~np.eye(N, dtype=bool)
        both = np.maximum(C, Cb)
        rmax = np.where(off, both, 0).max(1)
        rmin = float(np.minimum(C, Cb)[off].min())
        need = rmax > rmin + eps                     # one value for every pair: the curve alone is complete
        q_all = np.nonzero(need[L])[0]
        if not len(q_all):
            return []
        tree = cKDTree(P)
        cand = []
        for r in np.unique(rmax[L[q_all]]):
            q = q_all[rmax[L[q_all]] == r]
            r = min(float(r), self.max_gap)
            kk = int(min(len(P), 8 * r / self.step + 16))
            chunk = max(500, int(3e6 // kk))
            for s in range(0, len(q), chunk):
                qq = q[s:s + chunk]
                dist, idx = tree.query(P[qq], k=kk, distance_upper_bound=r, workers=-1)
                dist, idx = dist.reshape(len(qq), -1), idx.reshape(len(qq), -1)
                ok = idx < len(P)
                j = np.where(ok, idx, 0)
                ok &= L[j] != L[qq][:, None]
                if not ok.any():
                    continue
                ii, kj = np.nonzero(ok)
                p_i, j_i, d_i = qq[ii], j[ii, kj], dist[ii, kj]
                c = np.where(X[p_i] | X[j_i], Cb[L[p_i], L[j_i]], C[L[p_i], L[j_i]])
                v = d_i < c - eps
                p_i, j_i, d_i, c = p_i[v], j_i[v], d_i[v], c[v]
                if len(p_i):
                    _, f = np.unique(p_i.astype(np.int64) * (N + 2) + L[j_i], return_index=True)
                    cand.append(np.stack([p_i[f], j_i[f], d_i[f], c[f]], 1))
        if not cand:
            return []
        cand = np.vstack(cand)
        p_i, j_i, d_i, cc = cand[:, 0].astype(int), cand[:, 1].astype(int), cand[:, 2], cand[:, 3]
        # drop what the curve already shows for the same pair nearby (0.3 mm cells and their neighbours)
        la, lb, a, b = seen
        if len(la):
            lo, hi = np.minimum(la, lb), np.maximum(la, lb)
            keys = np.unique(_cell_key(np.concatenate([lo, lo]), np.concatenate([hi, hi]),
                                       P[np.concatenate([a, b])], N))
            lo, hi = np.minimum(L[p_i], L[j_i]), np.maximum(L[p_i], L[j_i])
            shown = np.zeros(len(p_i), bool)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    shown |= np.isin(_cell_key(lo, hi, P[p_i], N, dx, dy), keys)
            p_i, j_i, d_i, cc = p_i[~shown], j_i[~shown], d_i[~shown], cc[~shown]
        if not len(p_i):
            return []
        # a clear gap: shortened by eps at each end, it touches no copper
        u = (P[j_i] - P[p_i]) / np.maximum(d_i, 1e-12)[:, None]
        s0, s1 = P[p_i] + eps * u, P[j_i] - eps * u
        hit = np.zeros(len(p_i), bool)
        ci, _ = shapely.STRtree(self.geoms).query(shapely.linestrings(np.stack([s0, s1], 1)), predicate="intersects")
        hit[np.unique(ci)] = True
        p_i, j_i, d_i, cc = p_i[~hit], j_i[~hit], d_i[~hit], cc[~hit]
        if not len(p_i):
            return []
        # one violation per pair and neighbourhood
        mids = (P[p_i] + P[j_i]) / 2
        prs = np.minimum(L[p_i], L[j_i]) * (N + 2) + np.maximum(L[p_i], L[j_i])
        near = cKDTree(mids).query_pairs(LEVEL_RADIUS, output_type="ndarray")
        near = near[prs[near[:, 0]] == prs[near[:, 1]]] if len(near) else near
        n = len(p_i)
        gr = coo_matrix((np.ones(len(near)), (near[:, 0], near[:, 1])), shape=(n, n)) if len(near) \
            else coo_matrix((n, n))
        _, comp = connected_components(gr, directed=False)
        out = []
        for z in np.unique(comp):
            idx = np.nonzero(comp == z)[0]
            w = idx[np.argmax(cc[idx] / np.maximum(d_i[idx], 1e-12))]
            pts = np.vstack([P[p_i[idx]], P[j_i[idx]]])
            box = (pts[:, 0].min(), pts[:, 1].min(), pts[:, 0].max(), pts[:, 1].max())
            out.append(self._violation("level", int(L[p_i[w]]), int(L[j_i[w]]), d_i[w], cc[w], P[p_i[w]],
                                       P[j_i[w]], box, 0.0, X[p_i[w]] or X[j_i[w]]))
        return out

    # ------------------------------------------------------------------------------------------ records
    def _ordered(self, i, j, fa, fb):
        if self.nets[i] > self.nets[j]:
            i, j, fa, fb = j, i, fb, fa
        return i, j, (float(fa[0]), float(fa[1])), (float(fb[0]), float(fb[1]))

    def _violation(self, method, i, j, gap, required, fa, fb, box, length, bare):
        i, j, fa, fb = self._ordered(i, j, fa, fb)
        return Violation(self.nets[i], self.nets[j], float(gap), float(required), fa, fb, method, float(length),
                         tuple(float(v) for v in box), bool(bare))


def check(nets: Mapping[str, object], requirement: Requirement, step: float = 0.02, tau: float = 0.02,
          openings=None, bare_requirement: Requirement = None):
    """Build the curve of `nets` and check it once: see ClearanceCurve and ClearanceCurve.check."""
    return ClearanceCurve(nets, step, openings).check(requirement, bare_requirement, tau)


def _outline(geom):
    """The edges to sample: the boundary of an area, the geometry itself for lines and points."""
    if geom.geom_type in ("Point", "MultiPoint", "LineString", "MultiLineString", "LinearRing"):
        return geom
    return geom.boundary


def _cell_key(lo, hi, pts, N, dx=0, dy=0, cell=0.3):
    """One int64 per (net pair, 0.3 unit cell), the cell moved by dx, dy."""
    cx = np.floor(pts[:, 0] / cell).astype(np.int64) + dx + 2 ** 20
    cy = np.floor(pts[:, 1] / cell).astype(np.int64) + dy + 2 ** 20
    return ((lo.astype(np.int64) * (N + 2) + hi) * 2 ** 21 + cx) * 2 ** 21 + cy


def _zones(T, nb, band, viol):
    """Groups of violating ridges connected through shared triangles within the band (section 6)."""
    vi = np.nonzero(viol & band)[0]
    if not len(vi):
        return []
    bi = np.nonzero(band)[0]
    nodes, inv = np.unique(np.concatenate([T[bi], nb[bi]]), return_inverse=True)
    inv = inv.ravel()
    m = len(bi)
    gr = coo_matrix((np.ones(m), (inv[:m], inv[m:])), shape=(len(nodes), len(nodes)))
    _, comp = connected_components(gr, directed=False)
    rc = np.full(len(T), -1)
    rc[bi] = comp[inv[:m]]
    zc = rc[vi]
    order = np.argsort(zc, kind="stable")
    vi, zc = vi[order], zc[order]
    return np.split(vi, np.nonzero(np.diff(zc))[0] + 1)
