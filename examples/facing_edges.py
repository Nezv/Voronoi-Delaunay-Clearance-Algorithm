"""Two traces 0.15 mm apart, a pad near a plane, and a third net: one value for every pair."""
from shapely.geometry import Point, box

from clearance_curve import check

nets = {
    "A": box(0, 0, 5, 0.3),                     # two traces with a 0.15 mm gap
    "B": box(0, 0.45, 5, 0.75),
    "PAD": Point(2.5, 1.33).buffer(0.4),        # a round pad 0.18 mm above trace B
    "PLANE": box(0, 1.9, 5, 2.5),               # and a plane 0.17 mm above the pad
}

if __name__ == "__main__":
    for v in check(nets, requirement=0.2):
        print(f"{v.net_a:>5} - {v.net_b:<5}  gap {v.gap:.3f} mm  needs {v.required:.2f} mm  "
              f"over {v.length:.2f} mm of curve")
