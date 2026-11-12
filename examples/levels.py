"""Section 5: a 25 V net between a 35 V and a 0 V net hides their 0.54 mm neck from the curve.

The requirement here grows with the voltage between two nets: 0.6 mm above 30 V, 0.1 mm below (as IPC-2221B
Table 6-1, column B2, does for 31 to 100 V and 0 to 30 V). Use the figures of your own standard.
"""
from shapely.geometry import Point, box

from clearance_curve import ClearanceCurve

volts = {"35V": 35.0, "25V": 25.0, "0V": 0.0}
nets = {
    "35V": box(-3, 0.31, 3, 1).union(box(-0.1, 0.27, 0.1, 0.32)),          # a tab 0.27 mm above the middle
    "0V": box(-3, -1, 3, -0.31).union(box(-0.1, -0.32, 0.1, -0.27)),       # and one 0.27 mm below it
    "25V": box(-3, -0.15, -0.2, 0.15).union(Point(-0.2, 0).buffer(0.15)),  # a trace ending just before them
}


def requirement(a, b):
    return 0.6 if abs(volts[a] - volts[b]) > 30 else 0.1


if __name__ == "__main__":
    curve = ClearanceCurve(nets)
    plain = min(g.gap for g in curve.gaps(1.0) if {g.net_a, g.net_b} == {"35V", "0V"})
    print(f"narrowest 35V - 0V gap on the curve: {plain:.2f} mm")
    for v in curve.check(requirement):
        print(f"{v.net_a} - {v.net_b}: gap {v.gap:.2f} mm, needs {v.required:.2f} mm, found by the {v.method} search")
