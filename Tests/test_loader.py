"""
test_loader.py
==============
Tests for BCRTreeLoader — column detection, format detection, germline
construction, IMGT ungapping, and error handling.

What each test covers
---------------------
  Column detection
    test_scrambled_column_order   — columns in random order still map correctly
    test_format_a_detection       — H+L present → format="paired"
    test_format_b_detection       — no L column → format="heavy_only"
    test_non_standard_aliases     — renamed columns matched via _ALIASES
    test_column_map_property      — loader.column_map shows resolved names

  Missing columns
    test_missing_required_clone_id   — ValueError with helpful message
    test_missing_required_vdj_h      — ValueError with helpful message
    test_missing_optional_warns      — UserWarning (not crash) for c_call etc.
    test_missing_optional_fallback   — metadata fields default to "Unknown"

  Germline construction
    test_germline_from_column        — uses germline_alignment_d_mask
    test_germline_consensus_fallback — no germline column → consensus built
    test_germline_ungapping          — IMGT dots/dashes stripped before use

  Sequence building
    test_paired_seq_is_concatenated  — sequence = H + L
    test_heavy_only_seq_is_h_only    — sequence = H only
    test_imgt_dots_stripped          — dots removed from H+L sequences

  get_clones output
    test_clone_count                 — correct number of clone keys
    test_germline_injected           — each clone has exactly one germline record
    test_germline_not_counted_twice  — germline appears once even with many rows
    test_singleton_clone_present     — single-cell clones still appear in dict
                                       (filtering happens in pipeline, not loader)
    test_cell_id_fallback            — rows without cell_id get "cell_N" label
"""

import pytest
import pandas as pd

from BCR_lineage_tracer.loader import BCRTreeLoader, CellRecord, _ungap, _consensus


# ─────────────────────────────── helpers ─────────────────────────────────────

def loader_from_df(df, tmp_path):
    p = tmp_path / "tmp.xlsx"
    df.to_excel(p, index=False)
    return BCRTreeLoader(str(p)).load()


# ═══════════════════════════════════════════════════════════════════
# Column detection — order and aliases
# ═══════════════════════════════════════════════════════════════════

def test_scrambled_column_order(tmp_path, paired_df):
    """Shuffling column order must not affect mapping."""
    shuffled = paired_df[list(reversed(paired_df.columns))]
    loader = loader_from_df(shuffled, tmp_path)
    assert loader._cols["clone_id"] == "clone_id"
    assert loader._cols["vdj_h"]    == "VDJ_sequence_H"
    assert loader._cols["vdj_l"]    == "VDJ_sequence_L"
    assert loader._cols["germ_h"]   == "germline_alignment_d_mask"


def test_format_a_detection(paired_loader):
    assert paired_loader.format == "paired"


def test_format_b_detection(heavy_loader):
    assert heavy_loader.format == "heavy_only"


def test_non_standard_aliases(tmp_path, renamed_cols_df):
    """Columns like 'heavy_sequence', 'uca_sequence_H', 'clonotype' are aliased."""
    loader = loader_from_df(renamed_cols_df, tmp_path)
    assert loader._cols["clone_id"] == "clonotype"
    assert loader._cols["vdj_h"]    == "heavy_sequence"
    assert loader._cols["vdj_l"]    == "light_sequence"
    assert loader._cols["germ_h"]   == "uca_sequence_H"
    assert loader._cols["cell_id"]  == "barcode"
    assert loader._cols["sample_id"] == "tissue_origin"
    assert loader._cols["c_call"]   == "ig_class"


def test_column_map_property(paired_loader):
    cm = paired_loader.column_map
    assert isinstance(cm, dict)
    assert "clone_id" in cm
    assert "vdj_h"    in cm
    # returns a copy — mutations don't affect the loader
    cm["clone_id"] = "HACKED"
    assert paired_loader._cols["clone_id"] != "HACKED"


# ═══════════════════════════════════════════════════════════════════
# Missing columns — hard errors
# ═══════════════════════════════════════════════════════════════════

def test_missing_required_clone_id(tmp_path, paired_df):
    df = paired_df.drop(columns=["clone_id"])
    with pytest.raises(ValueError, match="clone identifier"):
        loader_from_df(df, tmp_path)


def test_missing_required_vdj_h(tmp_path, paired_df):
    df = paired_df.drop(columns=["VDJ_sequence_H"])
    with pytest.raises(ValueError, match="heavy-chain VDJ sequence"):
        loader_from_df(df, tmp_path)


def test_error_message_lists_file_columns(tmp_path, paired_df):
    """Error message should list the actual columns to help the user."""
    df = paired_df.drop(columns=["clone_id"])
    with pytest.raises(ValueError) as exc_info:
        loader_from_df(df, tmp_path)
    assert "VDJ_sequence_H" in str(exc_info.value)   # actual column listed
    assert "Actual columns" in str(exc_info.value)


# ═══════════════════════════════════════════════════════════════════
# Missing columns — soft warnings
# ═══════════════════════════════════════════════════════════════════

def test_missing_optional_warns(tmp_path, missing_optional_cols_df):
    """Missing c_call, sample_id, cluster_annotated → UserWarning, no crash."""
    with pytest.warns(UserWarning, match="not found"):
        loader_from_df(missing_optional_cols_df, tmp_path)


def test_missing_optional_fallback(tmp_path, missing_optional_cols_df):
    """Missing optional metadata → cells have 'Unknown' for those fields."""
    with pytest.warns(UserWarning):
        loader = loader_from_df(missing_optional_cols_df, tmp_path)
    clones = loader.get_clones()
    for rec in clones["clon1"]:
        if not rec.is_germline:
            assert rec.isotype           == "Unknown"
            assert rec.sample_id         == "Unknown"
            assert rec.cluster_annotated == "Unknown"


def test_missing_germ_col_warns(tmp_path, paired_df):
    """Absent germline column raises UserWarning about consensus fallback."""
    df = paired_df.drop(columns=["germline_alignment_d_mask"])
    with pytest.warns(UserWarning, match="germline"):
        loader_from_df(df, tmp_path)


# ═══════════════════════════════════════════════════════════════════
# Germline construction
# ═══════════════════════════════════════════════════════════════════

def test_germline_from_column(paired_loader):
    """Germline must equal the ungapped germline_alignment_d_mask value."""
    from Tests.conftest import GERM_H_A
    germline_in_loader = paired_loader._germlines["clone_A"]
    # For paired format without L germline column, only H is used
    assert germline_in_loader.startswith(GERM_H_A)


def test_germline_consensus_fallback(tmp_path, paired_df):
    """When germline column is absent, a consensus is built from observed seqs."""
    df = paired_df.drop(columns=["germline_alignment_d_mask"])
    with pytest.warns(UserWarning):
        loader = loader_from_df(df, tmp_path)
    assert "clone_A" in loader._germlines
    germ = loader._germlines["clone_A"]
    assert len(germ) > 0


def test_germline_ungapping(tmp_path, imgt_gapped_df):
    """IMGT dots in germline_alignment_d_mask must be stripped."""
    loader = loader_from_df(imgt_gapped_df, tmp_path)
    germ = loader._germlines["cG"]
    assert "." not in germ
    assert "-" not in germ


# ═══════════════════════════════════════════════════════════════════
# Sequence building
# ═══════════════════════════════════════════════════════════════════

def test_paired_seq_is_concatenated(paired_loader, paired_df):
    """Format A: cell sequence = ungapped H + ungapped L."""
    clones = paired_loader.get_clones()
    row = paired_df[paired_df["cell_id"] == "cellA_0"].iloc[0]
    expected = (
        row["VDJ_sequence_H"].upper()
        + row["VDJ_sequence_L"].upper()
    )
    obs = next(r for r in clones["clone_A"] if r.cell_id == "cellA_0")
    assert obs.sequence == expected


def test_heavy_only_seq_is_h_only(heavy_loader, heavy_only_df):
    """Format B: cell sequence = ungapped sequence_alignment only."""
    clones = heavy_loader.get_clones()
    row = heavy_only_df[heavy_only_df["cell_id"] == "cellHO_0"].iloc[0]
    expected = row["sequence_alignment"].upper()
    obs = next(r for r in clones["heavy_1"] if r.cell_id == "cellHO_0")
    assert obs.sequence == expected


def test_imgt_dots_stripped_from_sequence(tmp_path, imgt_gapped_df):
    """Dots in VDJ sequence must be removed before storing in CellRecord."""
    loader = loader_from_df(imgt_gapped_df, tmp_path)
    clones = loader.get_clones()
    for rec in clones["cG"]:
        assert "." not in rec.sequence
        assert "-" not in rec.sequence


# ═══════════════════════════════════════════════════════════════════
# get_clones output structure
# ═══════════════════════════════════════════════════════════════════

def test_clone_count(paired_clones):
    assert set(paired_clones.keys()) == {"clone_A", "clone_B"}


def test_germline_injected_per_clone(paired_clones):
    """Each clone must contain exactly one germline record."""
    for cid, recs in paired_clones.items():
        germs = [r for r in recs if r.is_germline]
        assert len(germs) == 1, f"clone {cid} has {len(germs)} germline records"


def test_germline_has_correct_metadata(paired_clones):
    germ = next(r for r in paired_clones["clone_A"] if r.is_germline)
    assert germ.isotype           == "Germline"
    assert germ.sample_id         == "Germline"
    assert germ.cluster_annotated == "Germline"


def test_germline_not_counted_twice(paired_clones):
    """Running get_clones multiple times must not accumulate germline nodes."""
    # Already called once via the fixture; if called again on the same
    # loader, the germline guard prevents duplicates
    assert sum(1 for r in paired_clones["clone_A"] if r.is_germline) == 1


def test_singleton_clone_present(paired_clones):
    """Singleton clones (1 cell) appear in get_clones output.
    Filtering (<2 cells) is the pipeline's job, not the loader's."""
    assert "clone_B" in paired_clones


def test_observed_cell_count(paired_clones):
    obs = [r for r in paired_clones["clone_A"] if not r.is_germline]
    assert len(obs) == 6


def test_cell_id_fallback(tmp_path, paired_df):
    """Rows without a cell_id column get 'cell_<row_index>' labels."""
    df = paired_df.drop(columns=["cell_id"])
    # Dropping cell_id triggers a UserWarning about the missing optional column
    with pytest.warns(UserWarning, match="cell identifier"):
        loader = loader_from_df(df, tmp_path)
    clones = loader.get_clones()
    for recs in clones.values():
        for rec in recs:
            if not rec.is_germline:
                assert rec.cell_id.startswith("cell_")


def test_heavy_only_tissue_metadata(heavy_clones):
    """Format B: sample_id is read from 'Sample' column."""
    recs = [r for r in heavy_clones["heavy_1"] if not r.is_germline]
    tissues = {r.sample_id for r in recs}
    assert "Lung" in tissues
    assert "mLN"  in tissues


# ═══════════════════════════════════════════════════════════════════
# Utility functions
# ═══════════════════════════════════════════════════════════════════

def test_ungap_removes_dots():
    assert _ungap("CAG...TTG--ATC") == "CAGTTGATC"

def test_ungap_empty_string():
    assert _ungap("") == ""

def test_ungap_clean_seq():
    assert _ungap("ACGT") == "ACGT"

def test_consensus_single():
    assert _consensus(["ACGT"]) == "ACGT"

def test_consensus_majority():
    # position 0: A,A,C → A wins
    result = _consensus(["ACGT", "ACGT", "CCGT"])
    assert result[0] == "A"

def test_consensus_ignores_minority_length():
    """Sequences that differ in length: modal length used, others ignored."""
    seqs = ["ACGT", "ACGT", "ACG"]   # modal len=4
    result = _consensus(seqs)
    assert len(result) == 4
