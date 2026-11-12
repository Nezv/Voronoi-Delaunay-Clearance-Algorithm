"""Section 7, the surface s: copper under solder mask may sit closer than bare copper. Gaps with a foot point in a
mask opening are checked against the bare requirement. gaps() measures every gap, violation or not."""
from shapely.geometry import box

from clearance_curve import ClearanceCurve

nets = {
    "VIN": box(0, 0, 6, 0.4),
    "SW": box(0, 0.6, 6, 1.0),                  # 0.2 mm from VIN along its whole length
}
openings = box(4, -0.2, 6.2, 1.2)               # the right third is bare (a pad's mask opening)

if __name__ == "__main__":
    curve = ClearanceCurve(nets, openings=openings)
    for v in curve.check(requirement=0.13, bare_requirement=0.25):
        print(f"{v.net_a} - {v.net_b}: gap {v.gap:.2f} mm, needs {v.required:.2f} mm, bare: {v.bare}, "
              f"x {v.box[0]:.2f} to {v.box[2]:.2f} mm")
    for g in curve.gaps(threshold=0.5):
        print(f"narrowest gap {g.net_a} - {g.net_b}: {g.gap:.2f} mm")
