"""
loader.py
=========
BCRTreeLoader — reads an Excel workbook and produces per-clone CellRecord
objects ready for tree construction.

Column detection
----------------
Column ORDER in the workbook does not matter.  Detection is name-based and
case-insensitive.  Each "role" (clone_id, vdj_h, …) has a list of aliases
tried in priority order; the first match found in the file wins.

Supported layouts
-----------------
  Format A  "paired"       e.g. PT7_data_for_course.xlsx
    Sequence  = VDJ_sequence_H  +  VDJ_sequence_L  (concatenated)
    Germline  = germline_alignment_d_mask (per row) or consensus fallback

  Format B  "heavy_only"   e.g. All_clones_heavy_isotype_clusters.xlsx
    Sequence  = sequence_alignment (heavy chain, IMGT-gapped → ungapped)
    Germline  = germline_alignment_d_mask

IMGT gaps
---------
Sequences stored in IMGT coordinate space contain '.' (framework insertions)
and '-' (deletions).  Both are stripped with _ungap() before any computation.
Germline sequences undergo the same treatment.

Missing columns
---------------
  clone_id, vdj_h → ValueError (cannot proceed without these)
  germ_h          → UserWarning; per-clone consensus germline built instead
  c_call,
  sample_id,
  cluster_annotated,
  cell_id         → graceful fallback to "Unknown" / row-index label;
                    a UserWarning lists every missing optional column so the
                    user knows what will be absent from the visualisation.
"""

from __future__ import annotations

import warnings
from collections import Counter
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .constants import GAP_CHARS


# ── helpers ──────────────────────────────────────────────────────────────────

def _first_col(columns: List[str], candidates: List[str]) -> Optional[str]:
    """Return the first candidate present in `columns` (case-insensitive)."""
    low = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand in columns:
            return cand
        if cand.lower() in low:
            return low[cand.lower()]
    return None


def _clean(val) -> str:
    if pd.isna(val):
        return ""
    return str(val).strip().upper()


def _ungap(seq: str) -> str:
    for ch in GAP_CHARS:
        seq = seq.replace(ch, "")
    return seq


def _consensus(seqs: List[str]) -> str:
    """Per-position modal base across a list of equal-length sequences."""
    if not seqs:
        return ""
    modal_len = Counter(len(s) for s in seqs).most_common(1)[0][0]
    modal_seqs = [s for s in seqs if len(s) == modal_len]
    arr = np.array([list(s) for s in modal_seqs])
    return "".join(
        Counter(arr[:, i]).most_common(1)[0][0] for i in range(arr.shape[1])
    )


# ── CellRecord ────────────────────────────────────────────────────────────────

class CellRecord:
    """One observed B cell (or synthetic germline) with its sequence and metadata.

    Attributes
    ----------
    cell_id           : unique identifier (barcode or row-index fallback)
    clone_id          : clonal family this cell belongs to
    sequence          : ungapped nucleotide sequence used for distance calc
    isotype           : heavy-chain constant-region call (c_call)
    sample_id         : tissue / sample of origin
    cluster_annotated : cell-type / GC-subset annotation
    is_germline       : True for the synthetic UCA / root node
    """

    __slots__ = ("cell_id", "clone_id", "sequence", "isotype",
                 "sample_id", "cluster_annotated", "is_germline")

    def __init__(self, cell_id: str, clone_id: str, sequence: str,
                 isotype: str = "Unknown", sample_id: str = "Unknown",
                 cluster_annotated: str = "Unknown", is_germline: bool = False):
        self.cell_id           = cell_id
        self.clone_id          = clone_id
        self.sequence          = sequence
        self.is_germline       = is_germline
        if is_germline:
            self.isotype           = "Germline"
            self.sample_id         = "Germline"
            self.cluster_annotated = "Germline"
        else:
            self.isotype           = isotype or "Unknown"
            self.sample_id         = sample_id or "Unknown"
            self.cluster_annotated = cluster_annotated or "Unknown"


# ── BCRTreeLoader ─────────────────────────────────────────────────────────────

class BCRTreeLoader:
    """Load a single-cell BCR xlsx table → per-clone CellRecord lists.

    Column order in the workbook is irrelevant.  Detection is name-based
    and case-insensitive via the _ALIASES table below.

    Usage
    -----
    >>> loader = BCRTreeLoader("data.xlsx").load()
    >>> print(loader.format)           # "paired" or "heavy_only"
    >>> print(loader.column_map)       # shows which file column maps to each role
    >>> clones = loader.get_clones()   # dict[clone_id → list[CellRecord]]
    """

    # ── alias table ──────────────────────────────────────────────────────────
    # Each key is an internal "role".  The list is tried left-to-right;
    # matching is case-insensitive.  Column ORDER in the file does not matter.
    #
    # To support a new naming convention, add the column name to the list.
    #
    # Required  : clone_id, vdj_h
    # Advised   : germ_h  (warning + consensus fallback if absent)
    # Optional  : everything else (graceful "Unknown" fallback)
    _ALIASES: Dict[str, List[str]] = {
        # identity
        "clone_id": [
            "clone_id", "CLONE_ID", "clone", "clonotype", "clonotype_id",
            "clone_group", "lineage_id",
        ],
        "cell_id": [
            "cell_id", "barcode", "sequence_id", "cell_barcode",
            "cell_name", "obs_name",
        ],
        # sequences
        "vdj_h": [
            "VDJ_sequence_H", "sequence_alignment_H", "sequence_alignment",
            "vdj_nt_H", "heavy_sequence", "nt_sequence_H",
        ],
        "vdj_l": [
            "VDJ_sequence_L", "sequence_alignment_L",
            "vdj_nt_L", "light_sequence", "nt_sequence_L",
        ],
        # germline  (IMGT-masked, D-region Ns)
        "germ_h": [
            "germline_alignment_d_mask", "germline_alignment_d_mask_H",
            "germline_alignment", "germline_sequence_H",
            "germline_nt_H", "uca_sequence_H",
        ],
        "germ_l": [
            "germline_alignment_d_mask_L",
            "germline_sequence_L", "germline_nt_L",
        ],
        # metadata
        "c_call": [
            "c_call", "C_CALL", "isotype", "c_gene", "constant_call",
            "heavy_isotype", "ig_class",
        ],
        "sample_id": [
            "sample_id", "Sample", "tissue", "organ", "origin",
            "tissue_origin", "tissue_of_origin", "biopsy_site",
        ],
        "cluster_annotated": [
            "cluster_annotated", "annot_clusters", "cell_type", "cluster",
            "cell_type_label", "annotation", "leiden", "seurat_clusters",
        ],
    }

    # Human-readable role descriptions used in error / warning messages
    _ROLE_DESC: Dict[str, str] = {
        "clone_id":          "clone identifier  (e.g. 'clone_id', 'clonotype')",
        "vdj_h":             "heavy-chain VDJ sequence  (e.g. 'VDJ_sequence_H', 'sequence_alignment')",
        "germ_h":            "IMGT germline  (e.g. 'germline_alignment_d_mask') — consensus fallback will be used",
        "c_call":            "isotype call  (e.g. 'c_call', 'isotype') — nodes will show 'Unknown'",
        "sample_id":         "tissue / sample  (e.g. 'sample_id', 'tissue') — nodes will show 'Unknown'",
        "cluster_annotated": "cell-type label  (e.g. 'cluster_annotated', 'cell_type') — nodes will show 'Unknown'",
        "cell_id":           "cell identifier  (e.g. 'cell_id', 'barcode') — row index used instead",
    }

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.df: Optional[pd.DataFrame] = None
        self.format: Optional[str] = None          # "paired" | "heavy_only"
        self._cols: Dict[str, Optional[str]] = {}  # role → actual column name
        self._germlines: Dict[str, str] = {}       # clone_id → germline seq

    # ── public API ────────────────────────────────────────────────────────────

    def load(self) -> "BCRTreeLoader":
        """Read workbook, detect format, pre-build germlines, emit warnings."""
        self.df = pd.read_excel(self.filepath)
        self._resolve_columns()
        self._detect_format()
        self._build_germlines()
        return self

    @property
    def column_map(self) -> Dict[str, Optional[str]]:
        """Show which actual file column was resolved for each role.

        Example output::

            {'clone_id': 'clone_id', 'vdj_h': 'sequence_alignment',
             'germ_h': 'germline_alignment_d_mask', 'c_call': 'c_call',
             'sample_id': 'Sample', 'cluster_annotated': 'annot_clusters', ...}
        """
        return dict(self._cols)

    def get_clones(self) -> Dict[str, List[CellRecord]]:
        """Return {clone_id: [CellRecord, …]} with one germline record per clone."""
        if self.df is None:
            raise RuntimeError("Call .load() before .get_clones()")

        clones: Dict[str, List[CellRecord]] = {}

        for idx, row in self.df.iterrows():
            raw = row.get(self._cols["clone_id"])
            if pd.isna(raw):
                continue
            cid = str(raw)

            seq = self._build_seq(row)
            if not seq:
                continue

            rec = CellRecord(
                cell_id=self._cell_id(row, idx),
                clone_id=cid,
                sequence=seq,
                isotype=self._get(row, "c_call") or "Unknown",
                sample_id=self._get(row, "sample_id") or "Unknown",
                cluster_annotated=self._get(row, "cluster_annotated") or "Unknown",
            )

            clones.setdefault(cid, [])

            # Inject one synthetic germline per clone
            gid = f"Germline_{cid}"
            if not any(r.cell_id == gid for r in clones[cid]):
                germ_seq = self._germlines.get(cid, "")
                if germ_seq:
                    clones[cid].append(CellRecord(
                        cell_id=gid, clone_id=cid, sequence=germ_seq,
                        is_germline=True,
                    ))

            clones[cid].append(rec)

        return clones

    # ── internal: column resolution ───────────────────────────────────────────

    def _resolve_columns(self):
        cols = list(self.df.columns)
        self._cols = {role: _first_col(cols, aliases)
                      for role, aliases in self._ALIASES.items()}

        # ── hard errors ──────────────────────────────────────────────────────
        missing_required = [
            role for role in ("clone_id", "vdj_h")
            if self._cols[role] is None
        ]
        if missing_required:
            lines = [
                "Could not find the following required column(s) in the workbook:\n"
            ]
            for role in missing_required:
                tried = self._ALIASES[role]
                lines.append(
                    f"  • {self._ROLE_DESC[role]}\n"
                    f"    Tried: {tried}\n"
                )
            lines.append(
                f"\nActual columns in the file:\n  {cols}\n\n"
                "Tip: add the correct column name to BCRTreeLoader._ALIASES['"
                + missing_required[0] + "'] or rename the column in your file."
            )
            raise ValueError("".join(lines))

        # ── soft warnings for optional-but-useful columns ─────────────────────
        missing_optional = [
            role for role in ("germ_h", "c_call", "sample_id",
                              "cluster_annotated", "cell_id")
            if self._cols[role] is None
        ]
        if missing_optional:
            lines = ["Some columns were not found and will use fallback values:\n"]
            for role in missing_optional:
                lines.append(f"  • {self._ROLE_DESC[role]}\n")
            lines.append(
                f"\nActual columns in the file:\n  {cols}\n"
                "If your file uses a different name, add it to "
                "BCRTreeLoader._ALIASES."
            )
            warnings.warn("".join(lines), UserWarning, stacklevel=3)

    def _detect_format(self):
        col = self._cols["vdj_l"]
        if col and self.df[col].apply(_clean).str.len().gt(0).any():
            self.format = "paired"
        else:
            self.format = "heavy_only"

    # ── internal: germline construction ───────────────────────────────────────

    def _build_germlines(self):
        """Populate self._germlines[clone_id] → ungapped germline sequence.

        Priority
        --------
        1. germline_alignment_d_mask column (IMGT-masked, per row)
           → modal-length ungapped sequences within each clone, consensus taken.
        2. Per-clone consensus of observed VDJ sequences (fallback when no
           germline column is present, e.g. PT7 file).
        """
        has_germ = self._cols["germ_h"] is not None
        col_cid  = self._cols["clone_id"]

        for cid, grp in self.df.groupby(col_cid):
            cid = str(cid)

            if has_germ:
                seqs = [_ungap(_clean(v))
                        for v in grp[self._cols["germ_h"]].dropna()]
                seqs = [s for s in seqs if s]
                if seqs:
                    self._germlines[cid] = _consensus(seqs)
                    continue

            # Fallback: consensus of observed sequences
            seqs_h = [_ungap(_clean(v))
                      for v in grp[self._cols["vdj_h"]].dropna()]
            seqs_h = [s for s in seqs_h if s]

            seqs_l: List[str] = []
            if self.format == "paired" and self._cols["vdj_l"]:
                seqs_l = [_ungap(_clean(v))
                          for v in grp[self._cols["vdj_l"]].dropna()]
                seqs_l = [s for s in seqs_l if s]

            if seqs_h:
                self._germlines[cid] = _consensus(seqs_h) + (
                    _consensus(seqs_l) if seqs_l else ""
                )

    # ── internal: per-row helpers ─────────────────────────────────────────────

    def _build_seq(self, row) -> str:
        h = _ungap(_clean(row.get(self._cols["vdj_h"], "")))
        if self.format == "paired" and self._cols["vdj_l"]:
            l = _ungap(_clean(row.get(self._cols["vdj_l"], "")))
            return h + l
        return h

    def _get(self, row, role: str) -> Optional[str]:
        col = self._cols.get(role)
        if not col:
            return None
        v = row.get(col)
        return None if pd.isna(v) else str(v)

    def _cell_id(self, row, idx) -> str:
        col = self._cols.get("cell_id")
        if col:
            v = row.get(col)
            if not pd.isna(v):
                return str(v)
        return f"cell_{idx}"
