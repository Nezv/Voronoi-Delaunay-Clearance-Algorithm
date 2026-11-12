"""Voronoi-Delaunay clearance: every gap between conductors closer than their requirement, found on the curve where
two conductors are equally near."""
from .curve import ClearanceCurve, Gap, Violation, check

__all__ = ["ClearanceCurve", "Gap", "Violation", "check"]
__version__ = "1.0.0"
