"""
app.py — FastAPI service for missense variant pathogenicity prediction.

Run:
    uvicorn app:app --reload --port 8000
Then open http://localhost:8000/docs  (interactive Swagger UI)

Design note (important ML-engineering principle): the service imports the SAME
features.py used in training. Training-serving skew — where features are computed
differently at train vs inference time — is one of the most common and pernicious
bugs in ML systems. Sharing one extraction function eliminates it by construction.
"""
from pathlib import Path
from typing import Optional

import joblib
import numpy as np
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

from features import extract_features

MODEL_PATH = Path(__file__).parent / "model.joblib"

app = FastAPI(
    title="Variant Effect Predictor",
    description="Predicts missense variant pathogenicity from biochemical "
                "sequence features.",
    version="0.1.0",
)

# Loaded once at startup, not per-request.
_bundle = None


def _get_model():
    global _bundle
    if _bundle is None:
        if not MODEL_PATH.exists():
            raise HTTPException(
                503,
                "model.joblib not found — train it first by running "
                "notebooks/03_modeling.ipynb.",
            )
        try:
            _bundle = joblib.load(MODEL_PATH)
        except Exception as e:  # corrupt file, sklearn version mismatch, etc.
            raise HTTPException(503, f"Failed to load model.joblib: {e}")
    return _bundle


class VariantRequest(BaseModel):
    variant: str = Field(..., examples=["p.Arg175His"],
                         description="HGVS protein change, 3- or 1-letter.")
    sequence: Optional[str] = Field(
        None, description="Optional full protein sequence for WT validation."
    )


class VariantResponse(BaseModel):
    variant: str
    score: float = Field(..., description="P(pathogenic), 0..1")
    prediction: str = Field(..., description="'pathogenic' or 'benign'")
    features: dict


@app.get("/health")
def health():
    return {"status": "ok", "model_loaded": MODEL_PATH.exists()}


@app.post("/predict", response_model=VariantResponse)
def predict(req: VariantRequest):
    bundle = _get_model()
    model, feature_cols = bundle["model"], bundle["feature_cols"]

    try:
        feats = extract_features(req.variant, sequence=req.sequence)
    except ValueError as e:
        # Bad variant string, WT mismatch, out-of-range position, etc.
        raise HTTPException(400, str(e))

    # Assemble the feature vector in the exact training column order.
    # Missing (sequence-dependent) features default to 0 — matches how the
    # model was trained (sequence-free), keeping train/serve consistent.
    x = np.array([[feats.get(c, 0.0) for c in feature_cols]])
    score = float(model.predict_proba(x)[0, 1])

    # Report ONLY the features the model actually consumed, so the response
    # never advertises sequence-derived features the sequence-free model
    # ignored. (extract_features returns 2 extra keys when a sequence is given.)
    model_feats = {c: feats.get(c, 0.0) for c in feature_cols}

    return VariantResponse(
        variant=req.variant,
        score=round(score, 4),
        prediction="pathogenic" if score >= 0.5 else "benign",
        features=model_feats,
    )