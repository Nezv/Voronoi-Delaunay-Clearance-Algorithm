"""Section 6: a wavy edge dips under the requirement eight times. With the band it is one violation."""
import numpy as np
from shapely.geometry import Polygon, box

from clearance_curve import ClearanceCurve

x = np.linspace(-1, 1, 2001)
edge = np.c_[x, 0.2 + 0.012 * np.sin(2 * np.pi * x / 0.25)]     # the gap swings from 0.188 to 0.212 mm
nets = {
    "A": box(-1, -0.5, 1, 0),
    "B": Polygon(np.vstack([edge, [(1, 0.6), (-1, 0.6)]])),
}

if __name__ == "__main__":
    curve = ClearanceCurve(nets, step=0.005)                    # built once, checked twice
    print("band 0.03 mm:", len(curve.check(0.2, tau=0.03)), "violation")
    print("no band:     ", len(curve.check(0.2, tau=0.0)), "violations")
