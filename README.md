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

All bioinformatics dependencies (BioPython, pandas, scikit-learn, XGBoost, and aligners) are defined in `environment.yml` and installed automatically during the build.

### Run the feature extractor

```bash
python features.py
```

The `extract_features` function parses HGVS protein notation, validates the wild-type residue against the reference sequence, and computes the full feature vector for each variant.

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

Labeled variants are sourced from **ClinVar** (variant summary) or the **dbNSFP** academic subset. Each record provides a pathogenic/benign label, HGVS protein notation, and an associated UniProt sequence. Mount your data into the container via the `./data` volume.

## Validation

Model performance is assessed via **ROC curves** with five-fold cross-validation.

## References

1. Adzhubei et al. (2010). *PolyPhen-2: Predicting the functional effect of human missense mutations.* Nature Methods.
2. Cheng et al. (2023). *Accurate proteome-wide missense variant effect prediction with AlphaMissense.* Science.

## License

Released under the MIT License.