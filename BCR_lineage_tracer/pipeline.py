"""
pipeline.py
===========
run()  —  top-level orchestration function that wires together the loader,
tracer, visualiser, and mutation-table exporter for every clone in a dataset.

Called by both the CLI (__main__.py) and the GUI (gui.py).

Parameters
----------
input_path          : path to the .xlsx workbook
output_dir          : directory to write tree PNGs and mutation_table.xlsx
clone_id            : if given, process only this one clone
collapse_threshold  : branch-length cutoff for polytomy collapsing
refine_isotypes     : enable isotype-aware NNI refinement
max_clones          : process at most this many clones (useful for debugging)
progress_cb         : optional callable(str) for live log messages
                      (GUI passes a function that appends to the log widget;
                       CLI leaves it None and we print() instead)

Returns
-------
(n_ok, n_skip, n_fail, combined_mutation_table)
"""

from __future__ import annotations

import itertools
import os
from typing import Callable, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from .loader import BCRTreeLoader, CellRecord
from .tracer import LineageTracer
from .visualization import plot_tree


def run(
    input_path: str,
    output_dir: str,
    clone_id: Optional[str] = None,
    collapse_threshold: float = 1e-6,
    refine_isotypes: bool = True,
    max_clones: Optional[int] = None,
    progress_cb: Optional[Callable[[str], None]] = None,
) -> Tuple[int, int, int, pd.DataFrame]:

    os.makedirs(output_dir, exist_ok=True)

    def log(msg: str):
        if progress_cb:
            progress_cb(msg)
        else:
            print(msg)

    # ── 1. Load data ────────────────────────────────────────────────────────
    loader = BCRTreeLoader(input_path).load()
    log(f"Format : {loader.format}")
    log(f"Rows   : {loader.df.shape[0]}")

    clones: Dict[str, List[CellRecord]] = loader.get_clones()
    log(f"Clones : {len(clones)} total")

    if clone_id:
        if clone_id not in clones:
            raise ValueError(f"clone_id '{clone_id}' not found in the input file.")
        clones = {clone_id: clones[clone_id]}

    if max_clones:
        clones = dict(itertools.islice(clones.items(), max_clones))

    # ── 2. Choose colour axis ───────────────────────────────────────────────
    # heavy-only → colour by tissue/organ (sample_id)
    # paired     → colour by isotype (c_call)
    color_by = "sample_id" if loader.format == "heavy_only" else "isotype"

    # ── 3. Process clones ───────────────────────────────────────────────────
    all_mutation_tables: List[pd.DataFrame] = []
    n_ok = n_skip = n_fail = 0

    for cid, records in clones.items():
        n_obs = sum(1 for r in records if not r.is_germline)

        if n_obs < 2:
            log(f"  SKIP  clone {cid}  ({n_obs} observed cell — need ≥ 2)")
            n_skip += 1
            continue

        try:
            tracer = LineageTracer(
                records,
                collapse_threshold=collapse_threshold,
                refine_isotypes=refine_isotypes,
            )
            tree = tracer.build()
        except Exception as exc:
            log(f"  FAIL  clone {cid}  — {exc}")
            n_fail += 1
            continue

        # ── visualisation
        out_png = os.path.join(output_dir, f"tree_{cid}.png")
        fig, _ = plot_tree(
            tree,
            color_by=color_by,
            shape_by="cluster_annotated",
            title=f"Clone {cid}  ({loader.format})",
            output_path=out_png,
        )
        plt.close("all")

        # ── mutation table
        mt = tracer.mutation_table()
        all_mutation_tables.append(mt)

        log(f"  ✓     clone {cid}  |  {n_obs} cells  |  {len(mt)} edges  →  {out_png}")
        n_ok += 1

    # ── 4. Export combined mutation table ───────────────────────────────────
    combined = (
        pd.concat(all_mutation_tables, ignore_index=True)
        if all_mutation_tables
        else pd.DataFrame()
    )
    mut_path = os.path.join(output_dir, "mutation_table.xlsx")
    combined.to_excel(mut_path, index=False)

    log("")
    log(f"Done  —  {n_ok} trees built,  {n_skip} skipped (<2 cells),  {n_fail} failed.")
    log(f"Mutation table → {mut_path}")

    return n_ok, n_skip, n_fail, combined
