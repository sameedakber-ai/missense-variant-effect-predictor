# Variant Effect Prediction

A machine-learning pipeline for predicting the pathogenicity of **missense variants** — single amino acid substitutions in human proteins — using chemical and evolutionary features.

## Overview

A missense variant swaps one amino acid for another. Some are harmless; others are catastrophic. This project quantifies *how disruptive a given substitution is at a given position* by combining two classical signals:

- **Chemical difference** between the original and substituted residue (charge, hydrophobicity, molecular weight)
- **Evolutionary acceptability** of the substitution, captured by the BLOSUM62 substitution matrix

These hand-crafted features feed a classifier that labels each variant as **pathogenic** or **benign** — the same features-plus-classifier architecture pioneered by PolyPhen-2.

## Background

The pipeline follows the central dogma: DNA → RNA → protein. Variants are differences from a reference genome; this project focuses exclusively on missense SNPs, where a base change alters a single residue in the resulting protein.

## Features

| Feature | Description |
|---------|-------------|
| Charge change | Difference in residue charge (positive/negative/neutral) |
| Hydrophobicity change | Shift in water affinity |
| Molecular weight change | Difference in residue size |
| BLOSUM62 score | Evolutionary substitution likelihood (the single most informative feature) |

## Project Structure

```text
variant_effect_prediction/
├─ Dockerfile              # Reproducible Miniconda-based environment
├─ docker-compose.yml      # Container orchestration and volume mounts
├─ environment.yml         # Conda/pip dependencies
├─ pyproject.toml          # Editable install so modules import in notebooks
├─ download_clinvar.py     # Fetch the ClinVar dataset with integrity check
├─ load_clinvar.py         # Load + filter + label raw ClinVar (memory-safe)
├─ clinvar_parse.py        # Parse ClinVar Names into protein consequences
├─ features.py             # Missense variant feature extraction
├─ app.py                  # FastAPI service exposing the trained model
├─ notebooks/              # Staged analysis (numbered, run in order)
│  ├─ 01_data_exploration.ipynb       # reads raw,       writes interim
│  ├─ 02_feature_engineering.ipynb    # reads interim,   writes processed
│  └─ 03_modeling.ipynb               # reads processed, writes model
├─ tests/                  # pytest for the module logic
├─ references/             # Source papers and notes
├─ data/                   # Datasets (git-ignored)
│  ├─ raw/                 # immutable downloaded source
│  ├─ interim/             # cleaned, labeled checkpoint
│  └─ processed/           # model-ready feature matrix
└─ results/                # Model outputs and figures
```

The pipeline is organized as three staged notebooks that hand off through
**Parquet checkpoints** on disk: each notebook reads one data layer and writes
the next, so any stage can be re-run independently. Pure logic (loading,
parsing, feature extraction) lives in importable modules; the notebooks
orchestrate and explore. The trained model is served over HTTP by `app.py`,
which imports the same `features.py` used in training so inference and training
compute features identically. See `README_ARCHITECTURE.md` for the full rationale.

## Getting Started

The project runs in a containerized, reproducible environment built with Docker and Miniconda. No local Python setup is required beyond Docker.

### Prerequisites

- [Docker](https://docs.docker.com/get-docker/) and Docker Compose

### Build and run

```bash
# Build the image and start the container
docker compose up -d --build

# Open a shell inside the running container
docker compose exec bioinfo bash
```

All bioinformatics dependencies (BioPython, pandas, scikit-learn, fastparquet, and aligners) are defined in `environment.yml` and installed automatically during the build.

### Enable clean notebook diffs (after cloning)

Notebooks are committed with their outputs stripped via [nbstripout](https://github.com/kynan/nbstripout), keeping the repository free of noisy, conflict-prone diffs. The filter is configured in the container automatically; if you also commit from the host, enable it there once:

```bash
nbstripout --install --attributes .gitattributes
```

### Download the data

```bash
python download_clinvar.py
```

This fetches the ClinVar `variant_summary` dataset into `./data/raw` and verifies the download against NCBI's published MD5 checksum. The file is gzipped and read directly — no manual decompression needed.

### Run the pipeline

The notebooks form an ordered pipeline; run them in sequence. Each reads the previous stage's checkpoint and writes its own:

```text
01_data_exploration.ipynb     raw        → data/interim/labeled_snvs.parquet
02_feature_engineering.ipynb  interim    → data/processed/feature_matrix.parquet
03_modeling.ipynb             processed  → model.joblib
```

`01` loads the ClinVar data via `load_clinvar.py` and reduces it to the confidently-labeled SNVs used downstream; `02` classifies consequences and builds the feature matrix; `03` trains and validates the classifier.

### Run the tests

```bash
pytest tests/
```

### Serve predictions (API)

Once `03_modeling.ipynb` has produced `model.joblib`, the prediction API starts
automatically with `docker compose up`. It runs alongside the `bioinfo`
container, reusing the same image so serving shares the training dependencies
and feature code.

```bash
# The api service is part of the compose stack; bring it up (or it's already up)
docker compose up -d api

# Interactive Swagger UI — try the model from the browser
# http://localhost:8000/docs

# Or call it directly
curl -X POST http://localhost:8000/predict \
    -H "Content-Type: application/json" \
    -d '{"variant": "p.Arg175His"}'

# Liveness / model-loaded check
curl http://localhost:8000/health
```

The `/docs` page is an auto-generated, interactive interface where users can
submit a variant and see the predicted pathogenicity score without writing any
client code.

### Launch Jupyter (optional)

```bash
docker compose run --rm --service-ports bioinfo \
    jupyter notebook --ip=0.0.0.0 --no-browser --allow-root
```

### Shut down

```bash
docker compose down
```

## Data

Labeled variants are sourced from **ClinVar** (variant summary) or the **dbNSFP** academic subset. Each record provides a pathogenic/benign label, HGVS protein notation, and an associated UniProt sequence.

Run `download_clinvar.py` to fetch the ClinVar variant summary into `./data/raw`, which is mounted into the container and excluded from version control. ClinVar publishes a new release on the first Thursday of each month; record the access date for reproducibility.

## Validation

Two models are trained and compared — an interpretable logistic-regression baseline and a gradient-boosting classifier. Performance is assessed via **ROC and precision-recall curves** with five-fold stratified cross-validation, plus a gene-held-out check to guard against leakage.

## References

1. Adzhubei et al. (2010). *PolyPhen-2: Predicting the functional effect of human missense mutations.* Nature Methods.
2. Cheng et al. (2023). *Accurate proteome-wide missense variant effect prediction with AlphaMissense.* Science.

## License

Released under the MIT License.