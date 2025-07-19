"""
View magnification methods for VTKPlot.

This module contains various algorithms for automatically focusing
the view on important regions of the data.
"""

from . import statistical
from . import variability
from . import field_range
from . import bounding_box
from . import mesh_refinement

__all__ = [
    "statistical",
    "variability", 
    "field_range",
    "bounding_box",
    "mesh_refinement"
]