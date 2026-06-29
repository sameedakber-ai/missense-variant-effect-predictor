"""
clinvar_parse.py
================

Turn the ClinVar `Name` column into something features.py can consume.

A ClinVar Name looks like:

    NM_000546.6(TP53):c.524G>A (p.Arg175His)

We want the protein-level part inside the parentheses: `p.Arg175His`. But the
real column is messy, and most of your data loss on Day 2 will come from NOT
handling these cases:

  - synonymous:        p.Leu125=        (no amino acid change -> not missense)
  - nonsense/stop:     p.Arg213Ter  or  p.Arg213*   (truncation -> not missense)
  - frameshift:        p.Lys382fs
  - unknown effect:    p.?  or  p.Met1?
  - extension/init:    p.Met1ext  / p.Met1Val (start-loss)
  - 3-letter only:     ClinVar almost always uses 3-letter (Arg, not R)

This module extracts ONLY clean single-residue missense changes
(wt != mut, both standard amino acids) and reports, by category, what it
skipped — so you can see your data loss instead of silently eating it.
"""

import re
import pandas as pd

# Capture the p.(...) consequence at the end of a ClinVar Name.
_P_PATTERN = re.compile(r"\(p\.([^)]+)\)")

# A clean 3-letter missense: WtAaPosMutAa, e.g. Arg175His
_MISSENSE_3 = re.compile(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

# Stop / synonymous / fs / unknown markers we want to classify (not crash on).
_STOP_TOKENS = ("Ter", "*")


def extract_protein_change(name: str):
    """
    Pull the raw p.-change string (without the 'p.') out of a ClinVar Name.
    Returns None if there's no protein consequence at all.
    """
    if not isinstance(name, str):
        return None
    m = _P_PATTERN.search(name)
    return m.group(1) if m else None


def classify_protein_change(p_change: str):
    """
    Classify a raw protein-change string (e.g. 'Arg175His', 'Leu125=',
    'Arg213Ter', 'Lys382fs', '?').

    Returns (category, parsed) where category is one of:
      'missense'    -> parsed = (wt3, pos, mut3) ; wt != mut, both standard
      'synonymous'  -> ends in '=' (no change)
      'nonsense'    -> introduces a stop
      'frameshift'  -> 'fs'
      'unknown'     -> '?', start-loss, or anything we won't model
    """
    if p_change is None:
        return ("unknown", None)

    s = p_change.strip()

    if s.endswith("="):
        return ("synonymous", None)
    if "fs" in s:
        return ("frameshift", None)
    if any(tok in s for tok in _STOP_TOKENS):
        return ("nonsense", None)
    if "?" in s or "ext" in s or "dup" in s or "del" in s or "ins" in s:
        return ("unknown", None)

    m = _MISSENSE_3.match(s)
    if not m:
        return ("unknown", None)

    wt3, pos, mut3 = m.group(1), int(m.group(2)), m.group(3)
    if wt3 == mut3:
        return ("synonymous", None)
    return ("missense", (wt3, pos, mut3))


def add_protein_change_columns(df: pd.DataFrame, name_col: str = "Name") -> pd.DataFrame:
    """
    Given a labeled ClinVar DataFrame, add:
       p_change      raw string inside p.(...)
       consequence   missense / synonymous / nonsense / frameshift / unknown
       p_hgvs        normalized 'p.Arg175His' for missense rows (else NaN)

    Then print a breakdown so you SEE your data funnel.
    """
    out = df.copy()
    out["p_change"] = out[name_col].map(extract_protein_change)

    cats, hgvs = [], []
    for pc in out["p_change"]:
        cat, parsed = classify_protein_change(pc)
        cats.append(cat)
        hgvs.append(f"p.{parsed[0]}{parsed[1]}{parsed[2]}" if cat == "missense" else None)

    out["consequence"] = cats
    out["p_hgvs"] = hgvs

    print("Consequence breakdown of labeled SNVs:")
    print(out["consequence"].value_counts(dropna=False).to_string())
    n_missense = (out["consequence"] == "missense").sum()
    print(f"\n-> {n_missense:,} clean missense rows will go to feature extraction.")
    return out


if __name__ == "__main__":
    # Realistic ClinVar Name strings to prove the classifier behaves.
    samples = [
        "NM_000546.6(TP53):c.524G>A (p.Arg175His)",      # missense
        "NM_000546.6(TP53):c.215C>G (p.Pro72Arg)",       # missense
        "NM_000546.6(TP53):c.375G>A (p.Thr125=)",        # synonymous
        "NM_000546.6(TP53):c.637C>T (p.Arg213Ter)",      # nonsense
        "NM_007294.4(BRCA1):c.5266dupC (p.Gln1756fs)",   # frameshift
        "NM_000546.6(TP53):c.1A>G (p.Met1?)",            # unknown/start
        "NM_000546.6(TP53):c.215C>G",                    # no protein part
        "single nucleotide variant with no name format", # junk
    ]
    for name in samples:
        pc = extract_protein_change(name)
        cat, parsed = classify_protein_change(pc)
        print(f"{cat:11s} | p.{pc if pc else '-':14s} | {name}")