"""
visualization.py
================
plot_tree() — renders a BCR clonal lineage tree as a rectangular cladogram
using matplotlib (headless-safe Agg backend).

Visual design decisions
-----------------------
- No axes frame: all four spines are hidden so no black box encloses the tree
  and long leaf labels at the right edge are never clipped by a border line.
- Legend placement: the axes right boundary is capped at 70 % of the figure
  width, and the legend is anchored at axes x = 1.02 (2 % gap).  The figure
  canvas is widened automatically based on the number of legend entries, so
  the legend never overlaps the tree area even for clones with many distinct
  isotypes or cell types.
- Node colour encodes the `color_by` attribute (isotype for paired data,
  sample_id / tissue for heavy-only data).
- Node shape encodes `cluster_annotated` (cell-type / GC-subset).
- The germline root node is drawn larger (160 pt²) with a square marker.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")   # must be set before pyplot import

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio.Phylo.BaseTree import Clade, Tree

from .constants import MARKER_CYCLE


# ── layout ────────────────────────────────────────────────────────────────────

def _layout(tree: Tree) -> Tuple[Dict[Clade, float], Dict[Clade, float]]:
    """Return x (cumulative branch length) and y (leaf slot) dicts."""
    x_pos: Dict[Clade, float] = {}

    def _assign_x(cl: Clade, x: float) -> None:
        x_pos[cl] = x
        for child in cl.clades:
            _assign_x(child, x + (child.branch_length or 0.0))

    _assign_x(tree.root, 0.0)

    y_pos: Dict[Clade, float] = {}
    for i, leaf in enumerate(tree.get_terminals()):
        y_pos[leaf] = float(i)

    def _assign_y(cl: Clade) -> float:
        if cl in y_pos:
            return y_pos[cl]
        child_ys = [_assign_y(c) for c in cl.clades]
        y = sum(child_ys) / len(child_ys)
        y_pos[cl] = y
        return y

    _assign_y(tree.root)
    return x_pos, y_pos


def _build_color_map(tree: Tree, color_by: str) -> Dict[str, object]:
    """Map each unique value of `color_by` to a matplotlib colour."""
    vals = sorted({str(getattr(c, color_by, "?")) for c in tree.find_clades()})
    cmap = plt.get_cmap("tab10")
    color_map: Dict[str, object] = {}
    for i, v in enumerate(vals):
        if v == "Germline":
            color_map[v] = "black"
        elif v == "Ancestral":
            color_map[v] = "lightgrey"
        else:
            color_map[v] = cmap(i % 10)
    return color_map


def _build_shape_map(tree: Tree, shape_by: str) -> Dict[str, str]:
    """Map each unique value of `shape_by` to a matplotlib marker string."""
    vals = sorted({str(getattr(c, shape_by, "?")) for c in tree.find_clades()})
    return {v: MARKER_CYCLE[i % len(MARKER_CYCLE)] for i, v in enumerate(vals)}


# ── main public function ──────────────────────────────────────────────────────

def plot_tree(
    tree: Tree,
    color_by: str = "isotype",
    shape_by: Optional[str] = "cluster_annotated",
    title: str = "",
    output_path: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    fig: Optional[plt.Figure] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """Draw a rectangular cladogram for `tree`.

    Parameters
    ----------
    tree         : annotated Bio.Phylo Tree from LineageTracer.build()
    color_by     : clade attribute used for node fill colour
    shape_by     : clade attribute used for node marker shape (None = circles)
    title        : figure title
    output_path  : if given, save PNG here (150 dpi, tight bounding box)
    ax / fig     : optional pre-existing axes to draw into

    Returns
    -------
    (fig, ax) tuple
    """
    x_pos, y_pos = _layout(tree)
    n_leaves = len(tree.get_terminals())

    # ── figure sizing ─────────────────────────────────────────────────────
    # We need enough horizontal space for:
    #   • the tree itself (drawn inside axes)
    #   • leaf-label text that extends past the axes right edge
    #   • the legend panel to the right of the axes
    #
    # Strategy:
    #   - Keep the tree axes occupying the LEFT 70 % of the figure width.
    #     This reserves 30 % for labels + legend, regardless of figure size.
    #   - Widen the figure canvas based on the number of legend entries so
    #     the legend always fits without overlapping.

    if ax is None:
        n_color = len({str(getattr(c, color_by, "?"))
                       for c in tree.find_clades()})
        n_shape = (len({str(getattr(c, shape_by, "?"))
                        for c in tree.find_clades()})
                   if shape_by else 0)
        n_legend_rows = n_color + n_shape

        # Base tree panel: 10 inches; legend panel: min 3.5 in, grows with entries
        legend_panel = max(3.5, n_legend_rows * 0.22)
        tree_panel   = 10.0
        fig_width    = tree_panel + legend_panel
        fig_height   = max(3.0, 0.35 * n_leaves)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Shrink the axes so it only occupies the left portion of the figure.
    # right=0.68 means the axes box ends at 68 % of figure width.
    # The remaining 32 % is shared by leaf label overflow and the legend.
    fig.subplots_adjust(left=0.02, right=0.68, top=0.96, bottom=0.04)

    # ── draw branches (rectangular cladogram) ────────────────────────────
    for cl in tree.find_clades():
        if not cl.clades:
            continue
        px = x_pos[cl]
        child_ys = [y_pos[c] for c in cl.clades]
        ax.plot([px, px], [min(child_ys), max(child_ys)],
                color="grey", lw=1, zorder=1)
        for child in cl.clades:
            ax.plot([px, x_pos[child]], [y_pos[child], y_pos[child]],
                    color="grey", lw=1, zorder=1)

    # ── colour / shape maps ───────────────────────────────────────────────
    color_map = _build_color_map(tree, color_by)
    shape_map = _build_shape_map(tree, shape_by) if shape_by else {}

    # ── draw nodes ────────────────────────────────────────────────────────
    for cl in tree.find_clades():
        x, y    = x_pos[cl], y_pos[cl]
        cv      = str(getattr(cl, color_by, "?"))
        is_germ = getattr(cl, "is_germline", False)

        marker = (
            shape_map.get(str(getattr(cl, shape_by, "?")), "o")
            if shape_by
            else ("s" if is_germ else "o")
        )

        ax.scatter(
            x, y,
            s=160 if is_germ else 70,
            marker=marker,
            facecolor=color_map.get(cv, "grey"),
            edgecolor="black",
            linewidth=0.8,
            zorder=3,
        )

        # Label all leaf nodes and the germline root
        if cl.is_terminal() or is_germ:
            ax.text(x, y, f"  {cl.name or ''}",
                    va="center", ha="left", fontsize=6, zorder=4)

    # ── legend ────────────────────────────────────────────────────────────
    # bbox_to_anchor uses FIGURE coordinates (bbox_transform=fig.transFigure)
    # so the legend position is independent of how far tree labels extend.
    # x=0.70 = just right of the 68 % axes boundary; y=0.97 = near top.
    handles: List[Line2D] = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=color_map[v],
               markeredgecolor="black", markersize=8,
               label=f"{color_by}: {v}")
        for v in color_map
    ]
    if shape_by:
        for v, mk in shape_map.items():
            handles.append(
                Line2D([0], [0], marker=mk, color="w",
                       markerfacecolor="lightgrey",
                       markeredgecolor="black", markersize=8,
                       label=f"{shape_by}: {v}")
            )

    ax.legend(
        handles=handles,
        bbox_to_anchor=(0.70, 0.97),       # figure-level coords
        bbox_transform=fig.transFigure,    # ← key: use figure not axes space
        loc="upper left",
        fontsize=6,
        frameon=True,
        framealpha=0.9,
        edgecolor="#cccccc",
        borderaxespad=0,
    )

    # ── axes cosmetics ────────────────────────────────────────────────────
    # Remove ALL four spines → no black frame around the tree.
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xlabel("Cumulative mutation distance from germline", labelpad=8)
    ax.set_yticks([])
    ax.tick_params(left=False)
    ax.set_title(title, fontsize=9, pad=10)

    # ── save ──────────────────────────────────────────────────────────────
    if output_path and fig:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)

    return fig, ax