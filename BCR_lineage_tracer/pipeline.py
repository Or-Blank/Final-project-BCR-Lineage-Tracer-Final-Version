"""
pipeline.py
===========
run() — top-level orchestration: loader → tracer → visualiser → Excel export.

Outputs per run
---------------
  tree_<clone_id>.png   Cladogram with seq1/seq2/… node labels and a
                        numeric x-axis showing mutation distance.

  mutation_table.xlsx   Per-edge mutation events with a seq_label column
                        (seq1, seq2, … for observed cells; anc1, anc2, …
                        for inferred internal nodes) so every row in the
                        table maps directly to a labelled node on the tree.
                        The original cell_id is kept alongside seq_label
                        for cross-referencing with the source data.
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

    def log(msg: str) -> None:
        if progress_cb:
            progress_cb(msg)
        else:
            print(msg)

    # ── 1. Load ───────────────────────────────────────────────────────────
    loader = BCRTreeLoader(input_path).load()
    log(f"Format : {loader.format}")
    log(f"Rows   : {loader.df.shape[0]}")

    clones: Dict[str, List[CellRecord]] = loader.get_clones()
    log(f"Clones : {len(clones)} total")

    if clone_id:
        if clone_id not in clones:
            raise ValueError(
                f"clone_id '{clone_id}' not found in the input file.")
        clones = {clone_id: clones[clone_id]}

    if max_clones:
        clones = dict(itertools.islice(clones.items(), max_clones))

    # ── 2. Colour axis ────────────────────────────────────────────────────
    color_by = "sample_id" if loader.format == "heavy_only" else "isotype"

    # ── 3. Process each clone ─────────────────────────────────────────────
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

        # ── visualisation — returns {cell_id → seq_label} mapping
        out_png = os.path.join(output_dir, f"tree_{cid}.png")
        fig, _, cell_id_map = plot_tree(
            tree,
            color_by=color_by,
            shape_by="cluster_annotated",
            title=f"Clone {cid}  ({loader.format})",
            output_path=out_png,
        )
        plt.close("all")

        # Build a complementary map for internal (ancestral) nodes.
        # The visualization module numbers them anc1, anc2, … in pre-order;
        # we replicate that here so the mutation table gets the same labels.
        anc_counter = 1
        anc_name_map: Dict[str, str] = {}
        for cl in tree.find_clades(order="preorder"):
            if cl.is_terminal():
                continue
            if getattr(cl, "is_germline", False):
                if cl.name:
                    anc_name_map[cl.name] = "Germline"
            elif cl.name and cl.name not in anc_name_map:
                anc_name_map[cl.name] = f"anc{anc_counter}"
                anc_counter += 1

        # Full node-name → short label (observed cells + ancestral nodes)
        full_label_map: Dict[str, str] = {**anc_name_map, **cell_id_map}

        # ── mutation table ────────────────────────────────────────────────
        mt = tracer.mutation_table()

        # Only add seq_label if the table has rows (avoids KeyError on
        # empty DataFrames which have no columns at all).
        if not mt.empty:
            # seq_label placed right after the "node" column
            node_pos = mt.columns.get_loc("node") + 1
            mt.insert(
                node_pos,
                "seq_label",
                mt["node"].map(full_label_map).fillna(""),
            )

        all_mutation_tables.append(mt)

        log(f"  ✓     clone {cid}  |  {n_obs} cells  "
            f"|  {len(cell_id_map)} seq labels  "
            f"|  {len(mt)} edges  →  {out_png}")
        n_ok += 1

    # ── 4. Export mutation table ──────────────────────────────────────────
    combined_mut = (
        pd.concat(all_mutation_tables, ignore_index=True)
        if all_mutation_tables
        else pd.DataFrame()
    )
    mut_path = os.path.join(output_dir, "mutation_table.xlsx")
    combined_mut.to_excel(mut_path, index=False)

    log("")
    log(f"Done  —  {n_ok} trees built,  "
        f"{n_skip} skipped (<2 cells),  {n_fail} failed.")
    log(f"Mutation table  → {mut_path}")

    return n_ok, n_skip, n_fail, combined_mut