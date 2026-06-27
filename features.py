"""
Project 1 - Variant Effect Prediction
Day 1 milestone: feature extraction for a single missense variant.

A missense variant swaps one amino acid (the wild-type, WT) for another (the 
mutant, MUT) at a given position in a protein. This module quantifies "how big 
a chemical insult is this swap" - the core intuition behind classical 
pathogenicity predictors like PolyPhen-2.

Run: python features.py
"""

import re
from typing import Optional, Tuple, cast, Dict
from Bio.Align import substitution_matrices

# ---------------------------------------------------------------------------
# Amino acid property tables (the 20 standard residues).
# These are the physicochemical scales the biology section referred to.
# ---------------------------------------------------------------------------

# Kyte-Doolittle hydrophobicity. Positive = hydrophobic (water-avoiding, 
# usually buried in the protein core; negative = hydrophilic (water-loving, 
# usually on the protein surface). A big change here can mean a residue that belongs 
# in the core is now exposed (or vice versa) - destabilizing.
KD_HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Net charge at physiological pH. Introducing or removing a charge - especially
# burying charge in the hydrophobic core - is frequently destabilizing.
CHARGE = {
    'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0,
    'Q': 0, 'E': -1, 'G': 0, 'H': 0.5, 'I': 0,
    'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
    'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
}

# Molecular weight (Da). A large size change can cause steric clashes if the 
# residue is packed tightly, or leave a destabilizing cavity if it shrinks.
MOL_WEIGHT = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
    'Q': 146.15, 'E': 147.13, 'G': 75.07, 'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}

# 3-letter -> 1-letter, for parsing HGVS notation like p.Arg175His
THREE_TO_ONE = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
}

_BLOSUM62 = substitution_matrices.load("BLOSUM62")
_BLOSUM62_DICT = cast(Dict[Tuple[str, str], float], _BLOSUM62)


def parse_hgvs_protein(variant: str) -> tuple[str, int, str]:
    """
    Parse HGVS protein notation into (wt, position, mut) as 1-letter codes.

    Accepts both 3-letter ('p.Arg175His') and 1-letter ('p.R175H') forms,
    with or without the leading 'p.'.
    """
    v = variant.strip()
    if v.startswith("p."):
        v = v[2:]

    # 3-letter form, e.g. Arg175His
    m = re.fullmatch(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', v)
    if m:
        wt = THREE_TO_ONE[m.group(1)]
        pos = int(m.group(2))
        mut = THREE_TO_ONE[m.group(3)]
        return wt, pos, mut

    # 1-letter form, e.g. R175H
    m = re.fullmatch(r'([A-Z])(\d+)([A-Z])', v)
    if m:
        wt = m.group(1)
        pos = int(m.group(2))
        mut = m.group(3)
        return wt, pos, mut

    raise ValueError(f"Could not parse variant: {variant!r}")


def extract_features(variant: str, sequence: Optional[str] = None) -> dict[str, float]:
    """
    Compute biochemical-difference features for a missense variant.

    Parameters
    ----------
    variant: str
        HGVS-style protein change, e.g. 'p.Arg175His' or 'p.R175H'
    sequence: str, optional
        Full protein sequence (1-letter). If given, we sanity-check that the 
        WT residue actually sits at the stated position, and add a couple of 
        position-relative features. This is the kind of validation real
        pipelines do - a mismatch means your variant and your sequence 
        disagree, which is a silent, dangerous bug.

    Returns
    -------
    dict of feature_name -> float
    """
    wt, pos, mut = parse_hgvs_protein(variant)

    if sequence is not None:
        if pos < 1 or pos > len(sequence):
            raise ValueError(f"Position {pos} out of range for length {len(sequence)}")
        observed = sequence[pos - 1]
        if observed != wt:
            raise ValueError(
                f"WT mismatch: variant says {wt} at {pos}, sequence has {observed}"
            )

    feats = {
        # The single most informative classical feature: how evolutionarily 
        # tolerated is the exact swap? Negative => rarely seen => disruptive.
        "blosum62": float(_BLOSUM62_DICT[wt, mut]),
        # Signed and absolute physiochemical deltas (MUT minus WT).
        "d_hydrophobicity": KD_HYDROPHOBICITY[mut] - KD_HYDROPHOBICITY[wt],
        "abs_d_hydrophobicity": abs(KD_HYDROPHOBICITY[mut] - KD_HYDROPHOBICITY[wt]),
        "d_charge": CHARGE[mut] - CHARGE[wt],
        "abs_d_charge": abs(CHARGE[mut] - CHARGE[wt]),
        "d_weight": MOL_WEIGHT[mut] - MOL_WEIGHT[wt],
        "abs_d_weight": abs(MOL_WEIGHT[mut] - MOL_WEIGHT[wt]),
        # Did the swap flip a residue between charged and uncharged? (0/1)
        "charge_introduced_or_removed": float(
            (CHARGE[wt] == 0) != (CHARGE[mut] == 0)
        ),

    }

    if sequence is not None:
        feats["rel_position"] = pos / len(sequence)
        feats["seq_length"] = float(len(sequence))

    return feats


if __name__ == "__main__":
    # TP53 R175H — one of the most famous cancer hotspot mutations.
    # (TP53 is "the guardian of the genome"; this exact variant recurs in
    #  many human tumors, so it's a great known-pathogenic sanity check.)
    demo_seq = (
        "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
        "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLQSGTAK"
        "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
        "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
        "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
        "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
    )
    for variant in ["p.Arg175His", "R175H", "p.Pro72Arg"]:
        try:
            feats = extract_features(variant, sequence=demo_seq)
            print(f"\n{variant}")
            for k, v in feats.items():
                print(f"  {k:30s} {v:8.3f}")
        except ValueError as e:
            print(f"\n{variant}: {e}")
