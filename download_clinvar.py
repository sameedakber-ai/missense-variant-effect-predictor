"""
download_clinvar.py
===================

Utilities for downloading the ClinVar ``variant_summary`` dataset from NCBI to a
local directory, with built-in integrity verification.

Overview
--------
ClinVar is a freely available, public archive maintained by the U.S. National
Center for Biotechnology Information (NCBI) that aggregates reports of human
genetic variants and their clinically asserted relationships to disease. The
``variant_summary.txt.gz`` file is a gzipped, tab-delimited table containing
every variant in ClinVar that has a defined location on the genome, along with
selected metadata for each: the clinical significance classification
(e.g. Pathogenic, Likely pathogenic, Benign), HGVS sequence nomenclature, gene
symbol, variant type, molecular consequence, and genomic coordinates. For
machine-learning work on variant pathogenicity, this file is the primary source
of labeled training data — clinical significance provides the target labels and
the HGVS notation provides the variant identifiers from which features are
derived.

This module retrieves that file (and its companion checksum) and verifies that
the downloaded bytes match what NCBI published.

Data source
-----------
Files are fetched from NCBI's public distribution tree:

    https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

Although this location is historically referred to as an "FTP site," NCBI now
serves it over HTTPS. The HTTPS endpoint is used here because it is more
reliable, more firewall- and proxy-friendly, and better supported by modern
HTTP client libraries than the legacy FTP protocol.

Two files are downloaded:

    variant_summary.txt.gz       The gzipped, tab-delimited variant table.
    variant_summary.txt.gz.md5   A small text file containing the expected MD5
                                 hash of the data file, used for verification.

Release cadence
---------------
NCBI publishes an updated ``variant_summary.txt.gz`` on the first Thursday of
each month. The live filename is not versioned; prior monthly releases are
moved to an ``archive/`` subdirectory and renamed with a date stamp
(e.g. ``variant_summary_2018-11.txt.gz``). Because the live file is overwritten
in place each month, callers wanting a specific monthly snapshot should record
the download date or retrieve the corresponding archived file directly.

Integrity verification
-----------------------
The data file is large (hundreds of megabytes) and a truncated or corrupted
download can silently produce malformed data. To guard against this, the
companion ``.md5`` file is downloaded alongside the data file, and the MD5 hash
of the local data file is computed and compared against the published value. A
mismatch indicates a corrupted or incomplete download that should be retried.

Implementation notes
---------------------
- Downloads are streamed to disk in chunks rather than buffered fully in memory,
  keeping memory usage flat regardless of file size.
- HTTP error responses (e.g. 404, 500) are raised as exceptions rather than
  being silently written to disk as if they were valid data.
- Existing files are skipped by default to avoid redundant re-downloads; this
  can be overridden to force a fresh pull of the latest monthly release.
- MD5 hashing reads the file incrementally, so verification also remains
  memory-light.

Usage
-----
As a script::

    python download_clinvar.py

As an importable function::

    from download_clinvar import download_clinvar_files

    download_clinvar_files()                 # default destination, skip if present
    download_clinvar_files("./data")         # custom destination directory
    download_clinvar_files(force=True)        # force re-download of latest release

Notes
-----
- The default destination ``/data`` resides at the filesystem root and typically
  requires elevated permissions to create. If a ``PermissionError`` occurs,
  supply a writable path within the project directory instead (e.g. ``"./data"``).
- The downloaded file remains gzipped. It does not need to be manually
  decompressed for analysis; ``pandas.read_csv`` reads ``.gz`` files directly via
  its ``compression`` parameter.

Dependencies
------------
- requests (third-party): HTTP client used for streaming downloads.
- Standard library: hashlib, shutil, pathlib.
"""


import hashlib
import shutil
from pathlib import Path

import requests


# Base URL of ClinVar's "FTP" tree (served over HTTPS).
CLINVAR_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/"

# The files we want. The .md5 is the checksum that lets us verify integrity.
CLINVAR_FILES = [
    "variant_summary.txt.gz",
    "variant_summary.txt.gz.md5",
]


def download_clinvar_files(dest_dir: str = "/data", force: bool = False) -> list[Path]:
    """
    Download the ClinVar variant_summary files into a local directory.

    Parameters
    ----------
    dest_dir : str
        Folder to save into. Created if it doesn't exist. Defaults to "/data".
    force : bool
        If False (default), skip files that already exist locally.
        If True, re-download and overwrite them.

    Returns
    -------
    list[Path]
        Paths to the files that now exist on disk.
    """
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)   # make /data if it's missing

    saved_paths: list[Path] = []

    for filename in CLINVAR_FILES:
        url = CLINVAR_BASE + filename
        out_path = dest / filename

        # Skip work if we already have it and aren't forcing a refresh.
        if out_path.exists() and not force:
            print(f"[skip] {filename} already exists at {out_path}")
            saved_paths.append(out_path)
            continue

        print(f"[download] {url}")

        # stream=True so we don't load the whole ~300 MB file into RAM at once.
        with requests.get(url, stream=True, timeout=60) as response:
            response.raise_for_status()   # raise an error on HTTP 404/500/etc.

            # Write to disk in chunks as bytes arrive.
            with open(out_path, "wb") as f:
                shutil.copyfileobj(response.raw, f)

        size_mb = out_path.stat().st_size / (1024 * 1024)
        print(f"[done] saved {filename} ({size_mb:.1f} MB)")
        saved_paths.append(out_path)

    # If we got both the data file and its checksum, verify integrity.
    data_file = dest / "variant_summary.txt.gz"
    md5_file = dest / "variant_summary.txt.gz.md5"
    if data_file.exists() and md5_file.exists():
        _verify_md5(data_file, md5_file)

    return saved_paths


def _verify_md5(data_file: Path, md5_file: Path) -> bool:
    """
    Check that data_file's MD5 hash matches the expected value in md5_file.
    NCBI's .md5 file contains a line like:  <hash>  variant_summary.txt.gz
    """
    expected = md5_file.read_text().split()[0].strip()

    # Compute the MD5 of the downloaded file, reading in chunks to stay memory-light.
    h = hashlib.md5()
    with open(data_file, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    actual = h.hexdigest()

    if actual == expected:
        print(f"[verify] MD5 OK: {actual}")
        return True
    else:
        print(f"[verify] MD5 MISMATCH!\n  expected: {expected}\n  actual:   {actual}")
        return False


if __name__ == "__main__":
    download_clinvar_files()