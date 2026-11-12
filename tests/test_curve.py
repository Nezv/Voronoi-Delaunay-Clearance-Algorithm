"""python -m unittest discover tests  (from the repository root)"""
import math
import os
import sys
import unittest

import numpy as np
from shapely.geometry import LineString, Point, box

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "examples"))

from clearance_curve import ClearanceCurve, check  # noqa: E402


class LocalGap(unittest.TestCase):
    def test_facing_edges(self):
        v = check({"A": box(0, 0, 5, 0.3), "B": box(0, 0.45, 5, 0.75)}, 0.2)
        self.assertEqual(len(v), 1)
        self.assertAlmostEqual(v[0].gap, 0.15, places=6)
        self.assertAlmostEqual(v[0].ratio, 0.2 / 0.15, places=6)
        self.assertEqual((v[0].net_a, v[0].net_b, v[0].method), ("A", "B", "curve"))
        self.assertAlmostEqual(v[0].length, 5.0, delta=0.05)

    def test_passes_at_the_requirement(self):
        self.assertEqual(check({"A": box(0, 0, 5, 0.3), "B": box(0, 0.5, 5, 0.8)}, 0.2), [])

    def test_round_pad_pinch(self):
        """Section 2: the narrowest gap is where g = 2 rho, straight below the pad's centre."""
        v = check({"PLANE": box(-2, 0, 2, 1), "PAD": Point(0.3, -0.5).buffer(0.35, 256)}, 0.2)
        self.assertEqual(len(v), 1)
        self.assertAlmostEqual(v[0].gap, 0.15, delta=0.002)
        self.assertAlmostEqual(v[0].foot_b[0], 0.3, delta=0.02)        # PLANE foot, above the centre

    def test_lines_and_points(self):
        v = check({"A": LineString([(0, 0), (3, 0)]), "B": LineString([(0, 0.1), (3, 0.1)])}, 0.15)
        self.assertAlmostEqual(v[0].gap, 0.1, places=6)

    def test_one_net_or_nothing(self):
        self.assertEqual(check({"A": box(0, 0, 1, 1)}, 0.2), [])
        self.assertEqual(check({}, 0.2), [])


class Splits(unittest.TestCase):
    def test_three_nets(self):
        """Section 3: around a split each pair keeps its own gap."""
        r = 0.5
        nets = {n: Point(math.cos(t) * 1.2, math.sin(t) * 1.2).buffer(r, 256)
                for n, t in zip("ABC", (math.pi / 2, math.pi / 2 + 2 * math.pi / 3, math.pi / 2 + 4 * math.pi / 3))}
        gaps = ClearanceCurve(nets, step=0.01).gaps(2.0)
        self.assertEqual(sorted((g.net_a, g.net_b) for g in gaps), [("A", "B"), ("A", "C"), ("B", "C")])
        exact = 2.4 * math.sin(math.pi / 3) - 2 * r
        for g in gaps:
            self.assertAlmostEqual(g.gap, exact, delta=0.002)


class Levels(unittest.TestCase):
    def test_hidden_neck(self):
        """Section 5: the plain curve shows 0.62 mm between 35 V and 0 V; the level search finds the 0.54 mm neck."""
        import levels
        curve = ClearanceCurve(levels.nets)
        plain = min(g.gap for g in curve.gaps(1.0) if {g.net_a, g.net_b} == {"35V", "0V"})
        self.assertAlmostEqual(plain, 0.62, places=4)
        v = curve.check(levels.requirement)
        self.assertEqual(len(v), 1)
        self.assertEqual((v[0].net_a, v[0].net_b, v[0].method), ("0V", "35V", "level"))
        self.assertAlmostEqual(v[0].gap, 0.54, places=4)

    def test_no_search_for_one_value(self):
        import levels
        self.assertEqual(ClearanceCurve(levels.nets).check(0.1), [])

    def test_blocked_chord(self):
        """A pair behind other copper is not a gap through air."""
        nets = {"HV": box(0, 0.6, 3, 1), "GND": box(0, -1, 3, -0.6), "MID": box(-1, -0.4, 4, 0.4)}
        req = lambda a, b: 2.0 if {a, b} == {"HV", "GND"} else 0.1
        self.assertEqual(check(nets, req), [])


class Counting(unittest.TestCase):
    def test_band(self):
        """Section 6: eight dips, one violation with the band, eight without."""
        import wavy_edge
        curve = ClearanceCurve(wavy_edge.nets, step=0.005)
        self.assertEqual(len(curve.check(0.2, tau=0.03)), 1)
        self.assertEqual(len(curve.check(0.2, tau=0.0)), 8)
        self.assertAlmostEqual(curve.check(0.2, tau=0.03)[0].gap, 0.188, delta=0.002)


class Surface(unittest.TestCase):
    def test_bare_requirement(self):
        import openings
        curve = ClearanceCurve(openings.nets, openings=openings.openings)
        v = curve.check(0.13, 0.25)
        self.assertEqual(len(v), 1)
        self.assertTrue(v[0].bare)
        self.assertGreaterEqual(v[0].box[0], 3.9)
        self.assertEqual(curve.check(0.13), [])

    def test_ratio_matches_check(self):
        import facing_edges
        curve = ClearanceCurve(facing_edges.nets)
        seg, g = curve.edges()
        u = curve.ratio(0.2)
        self.assertEqual(seg.shape, (len(g), 2, 2))
        self.assertAlmostEqual(float(u.max()), max(v.ratio for v in curve.check(0.2)), places=6)

    def test_requirement_function_is_symmetric(self):
        calls = []
        req = lambda a, b: calls.append((a, b)) or 0.2
        check({"A": box(0, 0, 1, 1), "B": box(1.1, 0, 2, 1), "C": box(2.1, 0, 3, 1)}, req)
        self.assertEqual(len(calls), 3)


if __name__ == "__main__":
    unittest.main()
