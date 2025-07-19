"""
VTKPlot: A VTK-based visualization tool for Tecplot data files.

This package provides visualization capabilities for Tecplot data files
with multiple view magnification methods that simulate Tecplot 360's
automatic view focusing.
"""

__version__ = "1.0.0"
__author__ = "Dmitry Mikushin"

from .core import process_tecplot_file, parse_tecplot_layout

__all__ = ["process_tecplot_file", "parse_tecplot_layout"]