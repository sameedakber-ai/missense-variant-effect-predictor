# Variant Effect Prediction — Pipeline Architecture

A machine-learning pipeline for predicting the pathogenicity of **missense
variants** from chemical and evolutionary features. This document describes the
**notebook + module architecture** and how the stages hand off to one another.

## The design principle

Three jobs are kept strictly separate:

- **Modules** (`*.py`) hold the logic that must be *correct* — pure
  transformations that are imported, reused, and unit-tested.
- **Notebooks** (`notebooks/*.ipynb`) hold the *reasoning* — what the data looks
  like, what we concluded, the figures. Each notebook narrates one stage.
- **Disk checkpoints** (`data/interim`, `data/processed`) hold the *handoffs*
  between stages, so each notebook can be re-run independently without
  recomputing the ones before it.

A notebook earns its boundary when it produces a stable artifact the next stage
consumes. If a notebook's job can be stated as *"reads X, produces Y,"* it is a
real stage.

## Data flow

```
download_clinvar.py          ->  data/raw/variant_summary.txt.gz
                                          |
01_data_exploration.ipynb    reads raw   ->  data/interim/labeled_snvs.parquet
02_feature_engineering.ipynb reads interim ->  data/processed/feature_matrix.parquet
03_modeling.ipynb            reads processed ->  model.joblib
```

Each notebook depends only on the **checkpoint file** of the previous stage,
never on another notebook's in-memory state. Run them in order: 01 → 02 → 03.

## Layout

```text
variant_effect_prediction/
├─ data/
│  ├─ raw/         # immutable input — variant_summary.txt.gz (git-ignored)
│  ├─ interim/     # labeled_snvs.parquet      (notebook 01 writes)
│  └─ processed/   # feature_matrix.parquet    (notebook 02 writes)
├─ notebooks/
│  ├─ 01_data_exploration.ipynb     # reads raw,       writes interim
│  ├─ 02_feature_engineering.ipynb  # reads interim,   writes processed
│  └─ 03_modeling.ipynb             # reads processed, writes model.joblib
├─ load_clinvar.py    # assembly/SNV filter + label cleaning (extracted logic)
├─ clinvar_parse.py   # ClinVar Name -> protein consequence classification
├─ features.py        # missense -> biochemical feature vector
├─ download_clinvar.py
├─ tests/             # pytest for the pure logic in the modules above
└─ model.joblib       # trained classifier + feature_cols (notebook 03 writes)
```

## What lives where, and why

| Concern | Location | Reason |
|---------|----------|--------|
| Label-cleaning rules (GRCh38, SNV, pathogenic/benign) | `load_clinvar.py` | Reviewer-scrutinized, branchy → testable module |
| Protein-consequence classification | `clinvar_parse.py` | Pure logic, reused across notebooks |
| Feature computation | `features.py` | Imported by notebook **and** the FastAPI service |
| "What does the data look like?" | notebooks | Exploration / narration, not reused |
| Plots, class-balance checks, funnels | notebooks | Discovery, belongs in the lab record |

**Rule of thumb:** code that *discovers* something stays in the notebook; code
that *does* a deterministic transformation moves to a module the moment it is
(a) called more than once, (b) worth unit-testing, or (c) needed in production.

## Running

```bash
# 1. fetch the data
python download_clinvar.py            # -> data/raw/variant_summary.txt.gz

# 2. run the notebooks in order (Jupyter, or via jupytext/papermill)
#    01 -> interim, 02 -> processed, 03 -> model.joblib

# 3. tests for the module logic
pytest tests/
```

## Checkpoints use Parquet, not CSV

Intermediate files are Parquet: dtypes are preserved (no `dtype=str`
re-coercion on reload), files are smaller, and reads are far faster than CSV.
The raw NCBI download stays gzipped as distributed.

## Notebook ↔ script pairing (optional but recommended)

Each notebook has a paired `.py` (percent-format) representation via
[jupytext](https://jupytext.readthedocs.io/). The `.py` files diff cleanly in
git and are what you review; regenerate either side with
`jupytext --sync notebooks/<name>.ipynb`. Combined with `nbstripout` (already
configured in the container), this keeps the repository free of noisy notebook
output diffs.
```bash
jupytext --set-formats ipynb,py:percent notebooks/*.ipynb   # enable pairing once
```
