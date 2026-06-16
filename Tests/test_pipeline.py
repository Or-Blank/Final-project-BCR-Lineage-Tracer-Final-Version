"""
test_pipeline.py
================
Integration tests for the run() pipeline function.

What each test covers
---------------------
  test_pipeline_paired_runs         — end-to-end on paired data; returns counts
  test_pipeline_heavy_only_runs     — end-to-end on heavy-only data
  test_pipeline_creates_tree_png    — PNG files written to output_dir
  test_pipeline_creates_mutation_xlsx — mutation_table.xlsx created
  test_pipeline_skips_singletons    — 1-cell clones counted in n_skip
  test_pipeline_single_clone_filter — --clone-id processes exactly one clone
  test_pipeline_max_clones          — --max-clones limits processing
  test_pipeline_bad_clone_id_raises — unknown clone_id → ValueError
  test_mutation_table_has_data      — returned DataFrame is non-empty
  test_pipeline_no_isotype_refine   — flag respected without crash
"""

import os
import pytest
import pandas as pd

from BCR_lineage_tracer.pipeline import run


def test_pipeline_paired_runs(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    n_ok, n_skip, n_fail, _ = run(paired_xlsx, out)
    assert n_ok   == 1   # clone_A built
    assert n_skip == 1   # clone_B skipped (singleton)
    assert n_fail == 0


def test_pipeline_heavy_only_runs(heavy_only_xlsx, tmp_path):
    out = str(tmp_path / "out")
    n_ok, n_skip, n_fail, _ = run(heavy_only_xlsx, out)
    assert n_ok   == 1
    assert n_skip == 1
    assert n_fail == 0


def test_pipeline_creates_tree_png(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    run(paired_xlsx, out)
    pngs = [f for f in os.listdir(out) if f.endswith(".png")]
    assert len(pngs) == 1
    assert "clone_A" in pngs[0]


def test_pipeline_creates_mutation_xlsx(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    run(paired_xlsx, out)
    assert os.path.exists(os.path.join(out, "mutation_table.xlsx"))


def test_pipeline_skips_singletons(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    _, n_skip, _, _ = run(paired_xlsx, out)
    assert n_skip >= 1


def test_pipeline_single_clone_filter(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    n_ok, n_skip, n_fail, _ = run(paired_xlsx, out, clone_id="clone_A")
    assert n_ok == 1
    pngs = [f for f in os.listdir(out) if f.endswith(".png")]
    assert len(pngs) == 1


def test_pipeline_max_clones(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    n_ok, _, _, _ = run(paired_xlsx, out, max_clones=1)
    assert n_ok <= 1


def test_pipeline_bad_clone_id_raises(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    with pytest.raises(ValueError, match="not found"):
        run(paired_xlsx, out, clone_id="DOES_NOT_EXIST")


def test_mutation_table_has_data(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    _, _, _, df = run(paired_xlsx, out)
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_pipeline_no_isotype_refine(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    n_ok, _, n_fail, _ = run(paired_xlsx, out, refine_isotypes=False)
    assert n_ok == 1
    assert n_fail == 0


def test_pipeline_progress_callback(paired_xlsx, tmp_path):
    out = str(tmp_path / "out")
    messages = []
    run(paired_xlsx, out, progress_cb=messages.append)
    assert any("Format" in m for m in messages)
    assert any("Done" in m for m in messages)


def test_pipeline_output_dir_created(paired_xlsx, tmp_path):
    out = str(tmp_path / "brand_new_dir" / "nested")
    run(paired_xlsx, out)
    assert os.path.isdir(out)
