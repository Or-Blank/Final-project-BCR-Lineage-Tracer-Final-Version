"""
conftest.py — shared pytest fixtures for BCR Lineage Tracer tests.
All fixtures are in-memory / tmp_path so no patient data is needed.
"""
import random
import pandas as pd
import pytest

from BCR_lineage_tracer.loader import BCRTreeLoader, CellRecord
from BCR_lineage_tracer.tracer import LineageTracer

BASES = "ACGT"
random.seed(42)

def make_seq(n): return "".join(random.choices(BASES, k=n))
def mutate(seq, n):
    seq = list(seq)
    for pos in random.sample(range(len(seq)), min(n, len(seq))):
        seq[pos] = random.choice([b for b in BASES if b != seq[pos]])
    return "".join(seq)

GERM_H_A  = make_seq(120)
GERM_L_A  = make_seq(90)
GERM_H_B  = make_seq(100)
GERM_H_HO = make_seq(130)


@pytest.fixture()
def paired_df():
    """Paired H+L dataset (Format A). clone_A=6 cells, clone_B=1 cell (skip)."""
    isotypes = ["IGHM","IGHM","IGHG1","IGHG1","IGHA1","IGHG1"]
    clusters = ["GC_DZ","GC_LZ","MBC","PC","GC_DZ","MBC"]
    rows = [
        {"cell_id": f"cellA_{i}", "clone_id": "clone_A", "sample_id": "PT7",
         "cluster_annotated": clusters[i], "c_call": isotypes[i],
         "VDJ_sequence_H": mutate(GERM_H_A, i+1),
         "VDJ_sequence_L": mutate(GERM_L_A, i),
         "germline_alignment_d_mask": GERM_H_A}
        for i in range(6)
    ]
    rows.append({"cell_id":"cellB_0","clone_id":"clone_B","sample_id":"PT7",
                 "cluster_annotated":"Unknown","c_call":"IGHM",
                 "VDJ_sequence_H": mutate(GERM_H_B,2),
                 "VDJ_sequence_L": mutate(make_seq(90),2),
                 "germline_alignment_d_mask": GERM_H_B})
    return pd.DataFrame(rows)


@pytest.fixture()
def heavy_only_df():
    """Heavy-only multi-tissue dataset (Format B). heavy_1=5 cells, heavy_2=1 (skip)."""
    tissues  = ["Lung","mLN","Lung","BM","mLN"]
    isotypes = ["IGHM","IGHG2B","IGHG2C","IGHM","IGHG1"]
    clusters = ["B.GC","B.GC","B.GC","B.CD19CONTROL","B.GC"]
    rows = [
        {"cell_id": f"cellHO_{i}", "clone_id": "heavy_1",
         "Sample": tissues[i], "annot_clusters": clusters[i],
         "c_call": isotypes[i],
         "sequence_alignment": mutate(GERM_H_HO, i+2),
         "germline_alignment_d_mask": GERM_H_HO}
        for i in range(5)
    ]
    rows.append({"cell_id":"cellHO_5","clone_id":"heavy_2","Sample":"Lung",
                 "annot_clusters":"B.GC","c_call":"IGHM",
                 "sequence_alignment": mutate(GERM_H_HO,3),
                 "germline_alignment_d_mask": GERM_H_HO})
    return pd.DataFrame(rows)


@pytest.fixture()
def renamed_cols_df():
    """Same data but every column uses a non-standard alias name."""
    rows = [
        {"barcode": f"bc_{i}", "clonotype": "clon1", "tissue_origin": "spleen",
         "cell_type_label": "GC_DZ", "ig_class": "IGHM",
         "heavy_sequence": mutate(GERM_H_A, i+1),
         "light_sequence": mutate(GERM_L_A, i),
         "uca_sequence_H": GERM_H_A}
        for i in range(4)
    ]
    return pd.DataFrame(rows)


@pytest.fixture()
def missing_optional_cols_df():
    """Paired data with c_call, sample_id, cluster_annotated absent → 'Unknown'."""
    rows = [
        {"cell_id": f"cell_{i}", "clone_id": "clon1",
         "VDJ_sequence_H": mutate(GERM_H_A, i+1),
         "VDJ_sequence_L": mutate(GERM_L_A, i),
         "germline_alignment_d_mask": GERM_H_A}
        for i in range(4)
    ]
    return pd.DataFrame(rows)


@pytest.fixture()
def imgt_gapped_df():
    """Sequences with IMGT dots — loader must ungap them."""
    germ = "CAGGTT...CAGCTG" + GERM_H_A[:40] + "...TGAAG"
    seq  = "CAGGTT...CAGCTG" + mutate(GERM_H_A[:40], 3) + "...TGAAG"
    rows = [
        {"cell_id": f"c{i}", "clone_id": "cG",
         "VDJ_sequence_H": seq, "germline_alignment_d_mask": germ,
         "c_call": "IGHM", "sample_id": "spleen", "cluster_annotated": "GC"}
        for i in range(3)
    ]
    return pd.DataFrame(rows)


def _xlsx(df, tmp_path, name="data.xlsx"):
    p = tmp_path / name
    df.to_excel(p, index=False)
    return str(p)

@pytest.fixture()
def paired_xlsx(tmp_path, paired_df):       return _xlsx(paired_df, tmp_path, "paired.xlsx")
@pytest.fixture()
def heavy_only_xlsx(tmp_path, heavy_only_df): return _xlsx(heavy_only_df, tmp_path, "heavy.xlsx")
@pytest.fixture()
def renamed_cols_xlsx(tmp_path, renamed_cols_df): return _xlsx(renamed_cols_df, tmp_path, "renamed.xlsx")

def _loader(df, tmp_path):
    return BCRTreeLoader(_xlsx(df, tmp_path)).load()

@pytest.fixture()
def paired_loader(paired_df, tmp_path):     return _loader(paired_df, tmp_path)
@pytest.fixture()
def heavy_loader(heavy_only_df, tmp_path):  return _loader(heavy_only_df, tmp_path)
@pytest.fixture()
def paired_clones(paired_loader):           return paired_loader.get_clones()
@pytest.fixture()
def heavy_clones(heavy_loader):             return heavy_loader.get_clones()

@pytest.fixture()
def paired_tree(paired_clones):
    t = LineageTracer(paired_clones["clone_A"])
    return t, t.build()

@pytest.fixture()
def heavy_tree(heavy_clones):
    t = LineageTracer(heavy_clones["heavy_1"])
    return t, t.build()
