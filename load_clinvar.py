"""
load_clinvar.py
===============

Load the raw ClinVar ``variant_summary`` table and reduce it to confidently
labeled single-nucleotide variants (SNVs) on a single genome assembly.

Why this is a module (and not just notebook cells)
--------------------------------------------------
The filtering and label-cleaning rules below are exactly the kind of logic a
thesis reviewer will scrutinize: which assembly was kept, what counts as
"pathogenic," what was dropped. Logic with branches you could get subtly wrong
belongs where it can be unit-tested and reused, not buried in a notebook cell.
Notebook 01 imports ``load_and_label_clinvar`` and then *explores* its output;
the notebook narrates, this module decides.

The data funnel
---------------
Raw ClinVar contains the same variant under both GRCh37 and GRCh38, plus every
variant type (insertions, deletions, CNVs, ...) and every clinical significance
value (including "Uncertain significance," "Conflicting," "drug response," ...).
This module applies three filters, in order:

    1. Assembly  == GRCh38            (avoid double-counting the same variant)
    2. Type      == SNV               (a missense change is a single base swap)
    3. ClinicalSignificance -> {0, 1} (keep only confident benign / pathogenic)

Rows whose significance does not resolve to a confident binary label are
dropped. The function optionally returns the per-stage counts so a notebook can
display the funnel rather than silently eating the losses.
"""

from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Label vocabulary.
#
# Standard ClinVar cleaning: "Likely" calls are conventionally folded into their
# confident counterparts; the combined "Pathogenic/Likely pathogenic" review
# label is treated as pathogenic (and symmetrically for benign). Everything else
# -- "Uncertain significance," "Conflicting interpretations," "drug response,"
# "not provided," etc. -- is intentionally excluded, because it does not give a
# trustworthy training target.
# ---------------------------------------------------------------------------
PATHOGENIC_LABELS = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
}
BENIGN_LABELS = {
    "Benign",
    "Likely benign",
    "Benign/Likely benign",
}

ASSEMBLY = "GRCh38"
SNV_TYPE = "single nucleotide variant"


def significance_to_label(sig: str) -> float:
    """
    Map a ClinVar ClinicalSignificance string to a binary label.

    Returns
    -------
    1.0   if the call is (likely) pathogenic
    0.0   if the call is (likely) benign
    np.nan otherwise (uncertain, conflicting, drug-response, ...) -> to be dropped
    """
    if sig in PATHOGENIC_LABELS:
        return 1.0
    if sig in BENIGN_LABELS:
        return 0.0
    return np.nan


def load_and_label_clinvar(
    path: Union[str, Path],
    return_funnel: bool = False,
    chunksize: int = 200_000,
) -> Union[pd.DataFrame, tuple[pd.DataFrame, dict[str, int]]]:
    """
    Load raw ClinVar and reduce to confidently labeled GRCh38 SNVs.

    Parameters
    ----------
    path : str | Path
        Path to ``variant_summary.txt.gz`` (gzipped, tab-delimited).
    return_funnel : bool
        If True, also return a dict of per-stage row counts so the caller can
        display the data funnel.

    Returns
    -------
    pd.DataFrame
        Labeled SNVs with an added integer ``label`` column (0=benign,
        1=pathogenic). Every other original column is preserved.
    dict[str, int], optional
        Per-stage counts (only if ``return_funnel=True``).

    chunksize : int
        Number of rows read per chunk. The full file is never held in memory;
        each chunk is filtered to its surviving rows before the next is read.

    Notes
    -----
    Memory safety is the whole point of this implementation. The full ClinVar
    table is tens of millions of rows wide; reading it all at once -- especially
    as Python ``str`` objects -- can exhaust container memory and get the kernel
    OOM-killed (it dies with no traceback). Three things keep memory flat:

    1. ``usecols`` -- only the ~5 columns we actually need are read; the other
       ~25 columns are never materialized.
    2. Chunked reading -- the file is streamed in blocks of ``chunksize`` rows;
       each block is filtered down to the handful of surviving rows immediately,
       so peak memory is one chunk, not the whole file.
    3. Category/efficient dtypes are applied only to the small final frame.
    """
    path = Path(path)

    # Only these columns are needed downstream. Reading ~5 of ~30 columns is the
    # single biggest memory win and also speeds up parsing substantially.
    needed_cols = ["Type", "Name", "GeneSymbol", "ClinicalSignificance", "Assembly"]

    funnel: dict[str, int] = {"raw": 0, "grch38": 0, "snv": 0, "labeled": 0}
    kept_chunks: list[pd.DataFrame] = []

    reader = pd.read_csv(
        path,
        sep="\t",
        compression="infer",     # reads .gz directly
        usecols=needed_cols,     # <-- only materialize the columns we use
        dtype="string",          # pandas' efficient string dtype, not object
        chunksize=chunksize,     # <-- stream the file in blocks
    )

    for chunk in reader:
        funnel["raw"] += len(chunk)

        # 1. One genome assembly only.
        chunk = chunk[chunk["Assembly"] == ASSEMBLY]
        funnel["grch38"] += len(chunk)

        # 2. Single nucleotide variants only.
        chunk = chunk[chunk["Type"] == SNV_TYPE]
        funnel["snv"] += len(chunk)

        # 3. Confident binary label only. Map -> drop unlabeled in this chunk.
        label = chunk["ClinicalSignificance"].map(significance_to_label)
        chunk = chunk.assign(label=label).dropna(subset=["label"])
        funnel["labeled"] += len(chunk)

        if len(chunk):
            kept_chunks.append(chunk)

    # Concatenate only the small surviving rows from every chunk.
    if kept_chunks:
        df = pd.concat(kept_chunks, ignore_index=True)
        df["label"] = df["label"].astype(int)
    else:
        df = pd.DataFrame(columns=needed_cols + ["label"])

    if return_funnel:
        return df, funnel
    return df


def print_funnel(funnel: dict[str, int]) -> None:
    """Pretty-print the per-stage data funnel produced by load_and_label_clinvar."""
    print("ClinVar data funnel:")
    print(f"  raw rows                : {funnel['raw']:,}")
    print(f"  -> GRCh38               : {funnel['grch38']:,}")
    print(f"  -> single nucleotide    : {funnel['snv']:,}")
    print(f"  -> confidently labeled  : {funnel['labeled']:,}")
    if funnel["labeled"]:
        # Share of the SNV pool that survived label cleaning.
        kept = funnel["labeled"] / funnel["snv"] if funnel["snv"] else 0.0
        print(f"  ({kept:.1%} of GRCh38 SNVs had a confident benign/pathogenic call)")


if __name__ == "__main__":
    # Smoke test against the default container data path.
    default_path = Path(__file__).resolve().parent / "data" / "raw" / "variant_summary.txt.gz"
    if default_path.exists():
        df_labeled, funnel = load_and_label_clinvar(default_path, return_funnel=True)
        print_funnel(funnel)
        print(df_labeled["label"].value_counts().rename({0: "benign", 1: "pathogenic"}))
    else:
        print(f"No data file at {default_path}; run download_clinvar.py first.")