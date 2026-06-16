"""
BCR Lineage Tracer
==================
B cell receptor clonal lineage tree reconstruction from single-cell data.

Package layout
--------------
  bcr_lineage_tracer/
    __init__.py         – public API re-exports
    constants.py        – biological constants (isotype order, gap chars, etc.)
    loader.py           – BCRTreeLoader  (data loading, format detection, germline)
    tracer.py           – LineageTracer  (NJ tree, Fitch parsimony, isotype NNI)
    visualization.py    – plot_tree      (matplotlib cladogram)
    pipeline.py         – run()          (orchestrates loader → tracer → viz → export)
    gui.py              – launch_gui()   (tkinter front-end)
    __main__.py         – CLI entry-point
"""

from .constants import ISOTYPE_ORDER, GAP_CHARS, MARKER_CYCLE
from .loader import BCRTreeLoader, CellRecord
from .tracer import LineageTracer
from .visualization import plot_tree
from .pipeline import run

__all__ = [
    "ISOTYPE_ORDER", "GAP_CHARS", "MARKER_CYCLE",
    "CellRecord", "BCRTreeLoader",
    "LineageTracer",
    "plot_tree",
    "run",
]
