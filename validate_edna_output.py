#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["pandas"]
# ///
"""
OceanOmics eDNA Pipeline Output Validator
SOP: OcOm_B218

Usage:
    validate_edna_output.py <project_dir> [--partial] [--log LOG] [--sing2 PATH]

Flags:
    --partial   Check shared pipeline directories (default: all dirs)
    --log       Log file path (default: <project_dir>/validation_<timestamp>.log)
    --sing2     Path to singularity images dir (default: $SING2 env var)
"""

import argparse
import csv
import json
import logging
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path

import pandas as pd
import xml.etree.ElementTree as ET
import zipfile as _zipfile

PASS = "PASS"
WARN = "WARN"
FAIL = "FAIL"

# ── xlsx reader (raw XML — bypasses openpyxl dimension bug) ───────────────────

_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"

def _read_xlsx_sheet(xlsx_path: Path, sheet_name: str) -> list[list]:
    """
    Read all rows from a named sheet in an xlsx file.
    Returns a list of rows; each row is a list of cell values (str/None).
    Uses raw zipfile+XML to avoid openpyxl's read_only dimension bug on FAIRe files.
    """
    with _zipfile.ZipFile(xlsx_path) as z:
        # Build shared strings table
        with z.open("xl/sharedStrings.xml") as f:
            ss_tree = ET.parse(f)
        strings = [
            (si.find(f".//{{{_NS}}}t") or si.find(".//{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t")
             ).text or ""
            if si.find(f".//{{{_NS}}}t") is not None else ""
            for si in ss_tree.findall(f".//{{{_NS}}}si")
        ]

        # Find sheet index by name
        with z.open("xl/workbook.xml") as f:
            wb_tree = ET.parse(f)
        sheet_elem = next(
            (s for s in wb_tree.findall(f".//{{{_NS}}}sheet")
             if s.get("name") == sheet_name),
            None
        )
        if sheet_elem is None:
            raise KeyError(f"Sheet '{sheet_name}' not found in {xlsx_path.name}")
        sheet_id = sheet_elem.get("sheetId")

        with z.open(f"xl/worksheets/sheet{sheet_id}.xml") as f:
            ws_tree = ET.parse(f)

        def col_index(ref: str) -> int:
            """Convert column letter(s) to 0-based index, e.g. 'A'->0, 'B'->1."""
            letters = "".join(c for c in ref if c.isalpha())
            idx = 0
            for ch in letters:
                idx = idx * 26 + (ord(ch.upper()) - ord("A") + 1)
            return idx - 1

        def cell_value(c_elem):
            t = c_elem.get("t")
            v = c_elem.find(f"{{{_NS}}}v")
            if v is None:
                return None
            if t == "s":
                return strings[int(v.text)]
            return v.text

        rows_out = []
        for row_elem in ws_tree.findall(f".//{{{_NS}}}row"):
            cells = row_elem.findall(f"{{{_NS}}}c")
            if not cells:
                continue
            # Determine row width from last cell reference
            last_ref = cells[-1].get("r", "A1")
            width = col_index(last_ref) + 1
            row = [None] * width
            for c in cells:
                ref = c.get("r", "")
                if ref:
                    row[col_index(ref)] = cell_value(c)
            if any(v is not None for v in row):
                rows_out.append(row)

        return rows_out


def _xlsx_sheet_names(xlsx_path: Path) -> list[str]:
    """Return sheet names in order."""
    with _zipfile.ZipFile(xlsx_path) as z:
        with z.open("xl/workbook.xml") as f:
            tree = ET.parse(f)
    return [s.get("name") for s in tree.findall(f".//{{{_NS}}}sheet")]


def _xlsx_sheet_id_map(xlsx_path: Path) -> dict[str, str]:
    """Return {sheet_name: sheetId} for all sheets."""
    with _zipfile.ZipFile(xlsx_path) as z:
        with z.open("xl/workbook.xml") as f:
            tree = ET.parse(f)
    return {s.get("name"): s.get("sheetId")
            for s in tree.findall(f".//{{{_NS}}}sheet")}

# Tabs expected in every FAIRe xlsx
FAIRE_EXPECTED_TABS  = {"README", "projectMetadata", "sampleMetadata",
                         "experimentRunMetadata", "taxaRaw", "taxaFinal",
                         "otuRaw", "otuFinal"}
# Standard FAIRe format tabs — present but not flagged as unexpected
FAIRE_KNOWN_EXTRA_TABS = {"Drop-down values"}

# sam_data columns added internally by the pipeline that are not in FAIRe sampleMetadata
SAM_DATA_INTERNAL_COLS = {"discarded", "sample_type", "use_for_filter"}

# Control sample name substrings — expected to be in phyloseq but absent from FAIRe sampleMetadata
CONTROL_PATTERNS = ("_NTC_", "_ITC_", "_EB_", "_WC_", "_FC_")

PHYLOSEQ_SIF = "bioconductor-phyloseq:1.54.sif"
R_SCRIPT     = Path(__file__).parent / "check_phyloseq.R"


# ── logging ────────────────────────────────────────────────────────────────────

def setup_logging(log_path: Path) -> logging.Logger:
    log = logging.getLogger("edna_validator")
    log.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s  %(levelname)-5s  %(message)s",
                             datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    log.addHandler(fh)
    log.addHandler(ch)
    return log


# ── result tracking ────────────────────────────────────────────────────────────

class Results:
    def __init__(self, log: logging.Logger):
        self._items: list[tuple[str, str, str, str]] = []
        self._log = log

    def add(self, assay: str, directory: str, status: str, message: str) -> str:
        self._items.append((assay, directory, status, message))
        if status == PASS:
            self._log.debug(f"[{assay}] {directory}: {message}")
        elif status == WARN:
            self._log.warning(f"[{assay}] {directory}: {message}")
        else:
            self._log.error(f"[{assay}] {directory}: {message}")
        return status

    def for_assay(self, assay: str):
        return [i for i in self._items if i[0] == assay]

    def counts(self):
        fails  = [(a, d, m) for a, d, s, m in self._items if s == FAIL]
        warns  = [(a, d, m) for a, d, s, m in self._items if s == WARN]
        passes = [i for i in self._items if i[2] == PASS]
        return fails, warns, passes

    def all_items(self):
        return self._items


# ── file / dir helpers ─────────────────────────────────────────────────────────

def check_file(path: Path, results: Results, assay: str, directory: str) -> bool:
    if not path.exists():
        results.add(assay, directory, FAIL, f"MISSING: {path.name}")
        return False
    if path.stat().st_size == 0:
        results.add(assay, directory, FAIL, f"EMPTY: {path.name}")
        return False
    results.add(assay, directory, PASS, f"OK: {path.name}")
    return True


def check_dir(path: Path, results: Results, assay: str, directory: str,
              check_empty: bool = True) -> bool:
    if not path.exists() or not path.is_dir():
        results.add(assay, directory, FAIL, f"MISSING DIR: {path.name}/")
        return False
    if check_empty and not any(path.iterdir()):
        results.add(assay, directory, FAIL, f"EMPTY DIR: {path.name}/")
        return False
    results.add(assay, directory, PASS, f"OK DIR: {path.name}/")
    return True


def flag_unexpected(dir_path: Path, expected: set[str],
                    results: Results, assay: str, directory: str):
    if not dir_path.exists():
        return
    for item in dir_path.iterdir():
        if item.name not in expected:
            results.add(assay, directory, WARN, f"UNEXPECTED: {item.name}")


# ── expected file lists ────────────────────────────────────────────────────────

def shared_expected(p: str, a: str) -> dict[str, list[str]]:
    """Files that are always checked (shared with collaborators)."""
    return {
        "05-lca": [
            f"{p}_{a}_asv_nt_lca_with_fishbase_output.tsv",
            f"{p}_{a}_asv_nt_lulucurated_lca_with_fishbase_output.tsv",
            f"{p}_{a}_asv_nt_lulucurated_taxa_final.tsv",
            f"{p}_{a}_asv_nt_lulucurated_taxa_raw.tsv",
            f"{p}_{a}_asv_nt_taxa_final.tsv",
            f"{p}_{a}_asv_nt_taxa_raw.tsv",
            f"{p}_{a}_asv_curateddb_lca_with_fishbase_output.tsv",
            f"{p}_{a}_asv_curateddb_lulucurated_lca_with_fishbase_output.tsv",
            f"{p}_{a}_asv_curateddb_lulucurated_taxa_final.tsv",
            f"{p}_{a}_asv_curateddb_lulucurated_taxa_raw.tsv",
            f"{p}_{a}_asv_curateddb_taxa_final.tsv",
            f"{p}_{a}_asv_curateddb_taxa_raw.tsv",
        ],
        "06-aquamap": [
            f"{p}_{a}_aquamaps_asv_nt.csv",
            f"{p}_{a}_aquamaps_asv_curateddb.csv",
            f"{p}_{a}_aquamaps_asv_nt_lulucurated.csv",
            f"{p}_{a}_aquamaps_asv_curateddb_lulucurated.csv",
        ],
        "06-phyloseq": [
            f"{p}_{a}_asv_nt_final_taxa.tsv",
            f"{p}_{a}_asv_nt_flagged_phyloseq.rds",
            f"{p}_{a}_asv_nt_lulucurated_final_taxa.tsv",
            f"{p}_{a}_asv_nt_lulucurated_flagged_phyloseq.rds",
            f"{p}_{a}_asv_curateddb_final_taxa.tsv",
            f"{p}_{a}_asv_curateddb_flagged_phyloseq.rds",
            f"{p}_{a}_asv_curateddb_lulucurated_final_taxa.tsv",
            f"{p}_{a}_asv_curateddb_lulucurated_flagged_phyloseq.rds",
        ],
        "07-faire": [
            f"{p}_{a}_asv_nt_final_faire_metadata.xlsx",
            f"{p}_{a}_asv_nt_lulucurated_final_faire_metadata.xlsx",
            f"{p}_{a}_asv_curateddb_final_faire_metadata.xlsx",
            f"{p}_{a}_asv_curateddb_lulucurated_final_faire_metadata.xlsx",
        ],
        "07-multiqc": [
            f"{p}_{a}_multiqc_report.html",
        ],
        "07-proportional_filter": [
            f"{p}_{a}_asv_nt_faire_taxa_filtered.tsv",
            f"{p}_{a}_asv_nt_lulucurated_faire_taxa_filtered.tsv",
            f"{p}_{a}_asv_nt_lulucurated_OTU_filtered.tsv",
            f"{p}_{a}_asv_nt_lulucurated_phyloseq_filtered.rds",
            f"{p}_{a}_asv_nt_lulucurated_phyloseq_taxa_filtered.tsv",
            f"{p}_{a}_asv_nt_lulucurated_proportional_stats.txt",
            f"{p}_{a}_asv_nt_OTU_filtered.tsv",
            f"{p}_{a}_asv_nt_phyloseq_filtered.rds",
            f"{p}_{a}_asv_nt_phyloseq_taxa_filtered.tsv",
            f"{p}_{a}_asv_nt_proportional_stats.txt",
            f"{p}_{a}_asv_curateddb_faire_taxa_filtered.tsv",
            f"{p}_{a}_asv_curateddb_lulucurated_faire_taxa_filtered.tsv",
            f"{p}_{a}_asv_curateddb_lulucurated_OTU_filtered.tsv",
            f"{p}_{a}_asv_curateddb_lulucurated_phyloseq_filtered.rds",
            f"{p}_{a}_asv_curateddb_lulucurated_phyloseq_taxa_filtered.tsv",
            f"{p}_{a}_asv_curateddb_lulucurated_proportional_stats.txt",
            f"{p}_{a}_asv_curateddb_OTU_filtered.tsv",
            f"{p}_{a}_asv_curateddb_phyloseq_filtered.rds",
            f"{p}_{a}_asv_curateddb_phyloseq_taxa_filtered.tsv",
            f"{p}_{a}_asv_curateddb_proportional_stats.txt",
        ],
    }


def all_only_expected(p: str, a: str) -> dict[str, list[str]]:
    """Additional files only checked with --all flag."""
    return {
        "01-cutadapt": [],  # checked by directory presence only
        "01-fastqc":   [],  # checked by directory presence only
        "01-seqkit_stats": [
            "assigned_seqkit_stats.txt",
            "final_seqkit_stats.txt",
            f"{p}_{a}_prefilter_seqkit_stats.txt",
            "raw_seqkit_stats.txt",
            "unknown_seqkit_stats.txt",
        ],
        "02-dada2": [
            f"{p}_{a}_asv.fa",
            f"{p}_{a}_asv_final_table.tsv",
            f"{p}_{a}_asv_table.csv",
            f"{p}_{a}_lca_input.tsv",
            f"{p}_{a}_seq_tab.rds",
        ],
        "03-lulu": [
            f"{p}_{a}_asv_curated_table.tab",
            f"{p}_{a}_asv_lulu_map.tab",
            f"{p}_{a}_asv_match_list.txt",
            f"{p}_{a}_curated_asv.fa",
        ],
        "04-blast": [
            f"{p}_{a}_asv_nt_blastn_results.txt",
            f"{p}_{a}_asv_curateddb_blastn_results.txt",
            f"{p}_{a}_curated_asv_nt_blastn_results.txt",
            f"{p}_{a}_curated_asv_curateddb_blastn_results.txt",
        ],
        "04-ocomnbc": [
            f"{p}_{a}_asv_nt_lulucurated_ocom_nbc_output.tsv",
            f"{p}_{a}_asv_nt_ocom_nbc_output.tsv",
            f"{p}_{a}_asv_curateddb_lulucurated_ocom_nbc_output.tsv",
            f"{p}_{a}_asv_curateddb_ocom_nbc_output.tsv",
        ],
        "07-pipeline_info": [
            f"{p}_{a}_samplesheet.valid.csv",
            "software_versions.yml",
        ],
    }


# ── 05-lca ─────────────────────────────────────────────────────────────────────

LCA_FISHBASE_COLS = {"domain", "phylum", "class", "order", "family",
                     "genus", "species", "otu", "querycoverage", "%id"}

def check_lca_dir(lca_dir: Path, p: str, a: str, results: Results):
    d = "05-lca"
    exp = shared_expected(p, a)["05-lca"]

    for fname in exp:
        fpath = lca_dir / fname
        if not check_file(fpath, results, a, d):
            continue

        try:
            if "_taxa_final.tsv" in fname:
                df = pd.read_csv(fpath, sep="\t")
                if "seq_id" in df.columns:
                    dupes = df["seq_id"][df["seq_id"].duplicated()].unique().tolist()
                    if dupes:
                        results.add(a, d, FAIL,
                            f"Duplicate seq_id in {fname}: {dupes}")
                    else:
                        results.add(a, d, PASS, f"seq_id unique in {fname}")
                else:
                    results.add(a, d, WARN, f"seq_id column missing from {fname}")

            elif "_lca_with_fishbase_output.tsv" in fname:
                df = pd.read_csv(fpath, sep="\t", nrows=1)
                cols = set(df.columns.str.lower())
                missing = LCA_FISHBASE_COLS - cols
                if missing:
                    results.add(a, d, WARN,
                        f"Missing expected columns in {fname}: {sorted(missing)}")
                else:
                    results.add(a, d, PASS, f"Expected columns present in {fname}")
        except Exception as e:
            results.add(a, d, WARN, f"Could not parse {fname}: {e}")

    flag_unexpected(lca_dir, set(exp), results, a, d)


# ── 06-aquamap ─────────────────────────────────────────────────────────────────

def check_aquamap_dir(aq_dir: Path, p: str, a: str, results: Results):
    d = "06-aquamap"
    exp = set(shared_expected(p, a)["06-aquamap"])

    for fname in exp:
        fpath = aq_dir / fname
        if not check_file(fpath, results, a, d):
            continue

        if "no_aquamap" in fname:
            results.add(a, d, WARN,
                f"no_aquamap filename — check decimalLatitude/decimalLongitude: {fname}")

        try:
            df = pd.read_csv(fpath, index_col=0)

            dupes = df.index[df.index.duplicated()].unique().tolist()
            if dupes:
                results.add(a, d, FAIL,
                    f"Duplicate species in {fname}: {dupes}")
            else:
                results.add(a, d, PASS, f"No duplicate species in {fname}")

            if df.isna().all().all():
                results.add(a, d, WARN,
                    f"All values are NA in {fname} — check lat/lon data")
        except Exception as e:
            results.add(a, d, WARN, f"Could not parse {fname}: {e}")

    # Unexpected files (including no_aquamap variants not in expected list)
    if aq_dir.exists():
        for item in aq_dir.iterdir():
            if item.name not in exp:
                if "no_aquamap" in item.name:
                    results.add(a, d, WARN,
                        f"no_aquamap file — check lat/lon: {item.name}")
                else:
                    results.add(a, d, WARN, f"UNEXPECTED: {item.name}")


# ── samplesheet helper ────────────────────────────────────────────────────────

def _read_samplesheet_discarded(path: Path) -> dict[str, bool]:
    """
    Read the pipeline's valid samplesheet and return {samp_name: is_discarded}.

    The valid samplesheet only contains samples where discarded=False (discarded
    samples are excluded before the pipeline runs). So absence from this file
    means a sample was either explicitly discarded or excluded from the run.

    Sample names in the samplesheet have a trailing run-index suffix (e.g. _T1).
    These are stripped so names match the FAIRe/phyloseq convention.
    """
    if not path.exists():
        return {}
    result = {}
    try:
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                name = row.get("samp_name", "").strip()
                disc = row.get("discarded", "False").strip().lower() == "true"
                # Strip trailing _T<digits> run-index suffix
                name = re.sub(r"_T\d+$", "", name)
                if name:
                    result[name] = disc
    except Exception:
        pass
    return result


# ── 06-phyloseq ────────────────────────────────────────────────────────────────

PHYLOSEQ_TAXA_COLS = {"domain", "phylum", "class", "order", "family", "genus", "species", "lca"}

def check_phyloseq_dir(ps_dir: Path, p: str, a: str, faire_dir: Path,
                        results: Results, sing2: str):
    d = "06-phyloseq"
    exp = shared_expected(p, a)["06-phyloseq"]

    # Load the valid samplesheet once for this assay (used to classify
    # samples that appear in FAIRe but are absent from the phyloseq object)
    samplesheet_path = ps_dir.parent / "07-pipeline_info" / f"{p}_{a}_samplesheet.valid.csv"
    samplesheet_discarded = _read_samplesheet_discarded(samplesheet_path)

    for fname in exp:
        fpath = ps_dir / fname
        if not check_file(fpath, results, a, d):
            continue

        if "_final_taxa.tsv" in fname:
            try:
                df = pd.read_csv(fpath, sep="\t", nrows=1)
                cols = set(df.columns.str.lower())
                missing = PHYLOSEQ_TAXA_COLS - cols
                if missing:
                    results.add(a, d, WARN,
                        f"Missing expected columns in {fname}: {sorted(missing)}")
                else:
                    results.add(a, d, PASS, f"Expected columns present in {fname}")
            except Exception as e:
                results.add(a, d, WARN, f"Could not parse {fname}: {e}")

    # R-based checks for each .rds file
    if not sing2:
        results.add(a, d, WARN, "Skipping phyloseq .rds checks — SING2 not set")
    elif not R_SCRIPT.exists():
        results.add(a, d, WARN, f"Skipping phyloseq .rds checks — {R_SCRIPT} not found")
    else:
        sif = os.path.join(sing2, PHYLOSEQ_SIF)
        if not os.path.exists(sif):
            results.add(a, d, WARN, f"Skipping phyloseq .rds checks — container not found: {sif}")
        else:
            for fname in exp:
                if not fname.endswith("_flagged_phyloseq.rds"):
                    continue
                rds = ps_dir / fname
                if not rds.exists():
                    continue

                # Derive the matching FAIRe xlsx path
                faire_stem = fname.replace("_flagged_phyloseq.rds",
                                           "_final_faire_metadata.xlsx")
                faire_path = faire_dir / faire_stem

                r_checks = run_phyloseq_r(rds, faire_path if faire_path.exists() else None,
                                           sif, results, a, d)
                # r_checks returns (ps_samples, ps_cols) for downstream use
                if r_checks and faire_path.exists():
                    ps_samples, ps_cols = r_checks
                    compare_samples_with_faire(ps_samples, ps_cols, faire_path,
                                               fname, results, a, d,
                                               samplesheet_discarded)

    flag_unexpected(ps_dir, set(exp), results, a, d)


def run_phyloseq_r(rds: Path, faire_path, sif: str,
                   results: Results, assay: str, d: str):
    """Run check_phyloseq.R via singularity. Returns (sample_names, col_names) or None."""
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tf:
        out_json = tf.name

    # Bind all relevant dirs so singularity can access them
    bind_dirs = {str(rds.parent)}
    if faire_path:
        bind_dirs.add(str(faire_path.parent))
    bind_dirs.add(str(R_SCRIPT.parent))
    bind_str = ",".join(f"{d}:{d}" for d in bind_dirs)

    cmd = [
        "singularity", "exec",
        "--bind", bind_str,
        sif,
        "Rscript", str(R_SCRIPT),
        str(rds), out_json,
    ]

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if proc.returncode != 0:
            results.add(assay, d, FAIL,
                f"[{rds.name}] R check failed: {proc.stderr.strip()[:300]}")
            return None

        with open(out_json) as f:
            data = json.load(f)

        ps_samples = None
        ps_cols    = None

        for item in data:
            status  = item["status"]
            message = item["message"]
            if status == "SAMPLE_NAMES":
                ps_samples = message.split(",") if message else []
            elif status == "SAM_DATA_COLS":
                ps_cols = message.split(",") if message else []
            else:
                results.add(assay, d, status, f"[{rds.name}] {message}")

        return (ps_samples, ps_cols)

    except subprocess.TimeoutExpired:
        results.add(assay, d, WARN, f"[{rds.name}] R check timed out")
        return None
    except Exception as e:
        results.add(assay, d, WARN, f"[{rds.name}] R check error: {e}")
        return None
    finally:
        try:
            os.unlink(out_json)
        except OSError:
            pass


def compare_samples_with_faire(ps_samples: list, ps_cols: list, faire_path: Path,
                                 rds_name: str, results: Results, assay: str, d: str,
                                 samplesheet_discarded: dict | None = None):
    """Compare phyloseq sample names and columns against FAIRe sampleMetadata.

    samplesheet_discarded: {samp_name: is_discarded} from the valid samplesheet.
    The valid samplesheet only contains samples that entered the pipeline (i.e.
    discarded=False); absence from it means the sample was excluded/discarded.
    """
    if ps_samples is None:
        return
    if samplesheet_discarded is None:
        samplesheet_discarded = {}

    _CTRL_CATS = {"negative control", "positive control", "pcr standard"}

    try:
        sheet_names = set(_xlsx_sheet_names(faire_path))
        if "sampleMetadata" not in sheet_names:
            results.add(assay, d, WARN,
                f"[{rds_name}] sampleMetadata tab missing from FAIRe xlsx")
            return

        rows = _read_xlsx_sheet(faire_path, "sampleMetadata")

        # FAIRe format: rows starting with # are metadata; first non-# row is header
        header_row = next(
            (r for r in rows
             if r[0] is not None and not str(r[0]).startswith("#")),
            None
        )
        if header_row is None or str(header_row[0]).strip() != "samp_name":
            results.add(assay, d, WARN,
                f"[{rds_name}] Could not find samp_name header in FAIRe sampleMetadata")
            return

        header_idx  = rows.index(header_row)
        headers     = [str(v).strip() if v is not None else "" for v in header_row]
        header_lower = [h.lower() for h in headers]

        # Build samp_name → samp_category lookup from FAIRe
        samp_cat_col = header_lower.index("samp_category") if "samp_category" in header_lower else None
        faire_samp_cats: dict[str, str] = {}

        faire_samples: list[str] = []
        for r in rows[header_idx + 1:]:
            if r[0] is None or not str(r[0]).strip():
                continue
            name = str(r[0]).strip()
            faire_samples.append(name)
            if samp_cat_col is not None and samp_cat_col < len(r) and r[samp_cat_col] is not None:
                faire_samp_cats[name] = str(r[samp_cat_col]).strip().lower()

        ps_set     = set(ps_samples)
        faire_set  = set(faire_samples)
        only_ps    = ps_set - faire_set
        only_faire = faire_set - ps_set

        # ── samples in phyloseq but not in FAIRe ──────────────────────────────
        if only_ps:
            controls   = {s for s in only_ps if any(p in s for p in CONTROL_PATTERNS)}
            unexpected = only_ps - controls
            if controls:
                results.add(assay, d, PASS,
                    f"[{rds_name}] Control samples in phyloseq correctly absent from "
                    f"FAIRe sampleMetadata: {sorted(controls)}")
            if unexpected:
                results.add(assay, d, WARN,
                    f"[{rds_name}] Non-control samples in phyloseq but missing from "
                    f"FAIRe sampleMetadata: {sorted(unexpected)}")

        # ── samples in FAIRe but not in phyloseq ──────────────────────────────
        if only_faire:
            controls_ok:    set[str] = set()   # non-ITC controls — expected
            itc_warn:       set[str] = set()   # ITC absent — concerning
            discarded_ok:   set[str] = set()   # real samples excluded from pipeline
            genuinely_miss: set[str] = set()   # undiscarded real samples — unexpected

            for samp in only_faire:
                cat      = faire_samp_cats.get(samp, "")
                is_ctrl  = cat in _CTRL_CATS or any(p in samp for p in CONTROL_PATTERNS)
                is_itc   = "_ITC_" in samp

                if is_itc:
                    itc_warn.add(samp)
                elif is_ctrl:
                    controls_ok.add(samp)
                else:
                    # Real sample: check against valid samplesheet
                    if samp in samplesheet_discarded:
                        if samplesheet_discarded[samp]:
                            discarded_ok.add(samp)   # explicitly discarded
                        else:
                            genuinely_miss.add(samp) # in pipeline but missing from phyloseq
                    else:
                        # Not in valid samplesheet — likely excluded before pipeline ran
                        discarded_ok.add(samp)

            if controls_ok:
                results.add(assay, d, PASS,
                    f"[{rds_name}] Control/negative samples absent from phyloseq — "
                    f"expected (controls may have 0 reads): {sorted(controls_ok)}")
            if itc_warn:
                results.add(assay, d, WARN,
                    f"[{rds_name}] ITC positive control(s) absent from phyloseq — "
                    f"ITC should have reads, investigate: {sorted(itc_warn)}")
            if discarded_ok:
                results.add(assay, d, PASS,
                    f"[{rds_name}] Real sample(s) absent from phyloseq because they were "
                    f"excluded/discarded from the pipeline run (0 reads expected): "
                    f"{sorted(discarded_ok)}")
            if genuinely_miss:
                results.add(assay, d, WARN,
                    f"[{rds_name}] Real sample(s) in FAIRe but absent from phyloseq "
                    f"This could be because the fw_index and rv_index values are wrong. Or it could be that we need to resequence this data: "
                    f"{sorted(genuinely_miss)}")

        if not only_ps and not only_faire:
            results.add(assay, d, PASS,
                f"[{rds_name}] Sample names match between sam_data and FAIRe sampleMetadata")

        # Column name comparison (case-insensitive), excluding known internal columns
        if ps_cols:
            ps_cols_lower    = {c.lower() for c in ps_cols} - SAM_DATA_INTERNAL_COLS
            faire_cols_lower = {h.lower() for h in headers if h and h != "samp_name"}
            missing = faire_cols_lower - ps_cols_lower
            extra   = ps_cols_lower - faire_cols_lower
            if missing:
                results.add(assay, d, WARN,
                    f"[{rds_name}] FAIRe sampleMetadata columns absent from sam_data: "
                    f"{sorted(missing)}")
            if extra:
                results.add(assay, d, WARN,
                    f"[{rds_name}] sam_data columns absent from FAIRe sampleMetadata: "
                    f"{sorted(extra)}")
            if not missing and not extra:
                results.add(assay, d, PASS,
                    f"[{rds_name}] sam_data columns match FAIRe sampleMetadata columns")

    except Exception as e:
        results.add(assay, d, WARN,
            f"[{rds_name}] Could not compare samples with FAIRe xlsx: {e}")


# ── 07-faire ───────────────────────────────────────────────────────────────────

def check_faire_dir(faire_dir: Path, p: str, a: str, results: Results):
    d = "07-faire"
    exp = shared_expected(p, a)["07-faire"]

    for fname in exp:
        fpath = faire_dir / fname
        if not check_file(fpath, results, a, d):
            continue
        _check_faire_xlsx(fpath, p, a, d, results)

    flag_unexpected(faire_dir, set(exp), results, a, d)


def _check_faire_xlsx(fpath: Path, p: str, a: str, d: str, results: Results):
    fname = fpath.name
    try:
        sheet_id_map = _xlsx_sheet_id_map(fpath)
        sheet_names  = set(sheet_id_map.keys())

        # Expected tabs
        missing_tabs = FAIRE_EXPECTED_TABS - sheet_names
        if missing_tabs:
            results.add(a, d, FAIL,
                f"Missing tabs in {fname}: {sorted(missing_tabs)}")
        else:
            results.add(a, d, PASS, f"All expected tabs present in {fname}")

        # Unexpected tabs
        extra = sheet_names - FAIRE_EXPECTED_TABS - FAIRE_KNOWN_EXTRA_TABS
        if extra:
            results.add(a, d, WARN, f"Unexpected tabs in {fname}: {sorted(extra)}")

        # Project ID in projectMetadata
        # projectMetadata is transposed: col[2] = term_name, col[3] = value
        if "projectMetadata" in sheet_names:
            try:
                pm_rows = _read_xlsx_sheet(fpath, "projectMetadata")
                found = any(
                    len(r) > 3 and r[3] is not None and p in str(r[3])
                    for r in pm_rows
                )
                if not found:
                    results.add(a, d, WARN,
                        f"Project ID '{p}' not found in projectMetadata of {fname}")
                else:
                    results.add(a, d, PASS,
                        f"Project ID '{p}' found in projectMetadata of {fname}")
            except Exception as e:
                results.add(a, d, WARN, f"Could not read projectMetadata of {fname}: {e}")

        # Per-tab checks
        for tab in FAIRE_EXPECTED_TABS:
            if tab not in sheet_names:
                continue
            try:
                rows = _read_xlsx_sheet(fpath, tab)
            except Exception as e:
                results.add(a, d, WARN, f"Could not read tab '{tab}' of {fname}: {e}")
                continue

            if not rows:
                results.add(a, d, FAIL, f"Empty tab '{tab}' in {fname}")
                continue

            # Mandatory-column completeness check for all row-based tabs
            if tab in _FAIRE_ROW_TABS:
                sheet_id   = sheet_id_map.get(tab)
                conditions = _read_faire_conditions(fpath, sheet_id) if sheet_id else {}
                _check_mandatory_columns(rows, tab, fname, results, a, d, conditions)

            # taxaFinal: first column (ASV id) must be unique (skip metadata rows)
            if tab == "taxaFinal":
                data_rows = _skip_metadata_rows(rows)
                if data_rows:
                    first_col = [str(r[0]) for r in data_rows[1:] if r[0] is not None]
                    dupes = [k for k, v in Counter(first_col).items() if v > 1]
                    if dupes:
                        results.add(a, d, FAIL,
                            f"Duplicate ASV IDs in '{tab}' of {fname}: {dupes}")
                    else:
                        results.add(a, d, PASS,
                            f"ASV IDs unique in '{tab}' of {fname}")

            # otuRaw / otuFinal: sample columns (header row) must be unique
            # FAIRe OTU tab format: row 1 = ['', 'ASV', sample1, sample2, ...]
            # (no # metadata rows, so _skip_metadata_rows would land on data rows)
            if tab in ("otuRaw", "otuFinal"):
                if rows:
                    # First row is always the header; sample names start at col index 2
                    header = [str(v) for v in rows[0][2:] if v is not None]
                    dupes  = [k for k, v in Counter(header).items() if v > 1]
                    if dupes:
                        results.add(a, d, FAIL,
                            f"Duplicate sample columns in '{tab}' of {fname}: {dupes}")
                    else:
                        results.add(a, d, PASS,
                            f"Sample columns unique in '{tab}' of {fname}")

    except Exception as e:
        results.add(a, d, WARN, f"Could not inspect {fname}: {e}")


def _skip_metadata_rows(rows: list) -> list:
    """Skip FAIRe comment/metadata rows (those starting with # or requirement codes)."""
    for i, r in enumerate(rows):
        val = str(r[0]).strip() if r[0] else ""
        if not val.startswith("#") and val not in ("M", "R", "O", ""):
            return rows[i:]
    return rows


def _read_faire_conditions(xlsx_path: Path, sheet_id: str) -> dict[int, object]:
    """
    Parse the cell comments for a FAIRe sheet to extract conditional mandatory rules.

    FAIRe comment format (on each header-row cell):
        "Requirement level : Mandatory (If samp_category = negative control)"
        "Requirement level : Mandatory (Mandatory unless samp_category = negative control, positive control, or PCR standard)"

    Returns a dict  {col_index: condition_fn}  where condition_fn(samp_cat: str) -> bool
    returns True when the column IS mandatory for a row with that samp_category.
    Only columns with conditional (not absolute) M requirements are returned.
    Absolute M columns are handled by the regular mandatory-index list.
    """
    _NS_REL = "http://schemas.openxmlformats.org/package/2006/relationships"

    def _col_idx(ref: str) -> int:
        letters = "".join(c for c in ref if c.isalpha())
        idx = 0
        for ch in letters:
            idx = idx * 26 + (ord(ch.upper()) - ord("A") + 1)
        return idx - 1

    try:
        with _zipfile.ZipFile(xlsx_path) as z:
            rels_path = f"xl/worksheets/_rels/sheet{sheet_id}.xml.rels"
            if rels_path not in z.namelist():
                return {}
            with z.open(rels_path) as f:
                rels_tree = ET.parse(f)

            comments_file = None
            for rel in rels_tree.findall(f".//{{{_NS_REL}}}Relationship"):
                if rel.get("Type", "").endswith("/comments"):
                    target = rel.get("Target", "").lstrip("../")
                    comments_file = f"xl/{target}"
                    break

            if not comments_file or comments_file not in z.namelist():
                return {}

            with z.open(comments_file) as f:
                cm_tree = ET.parse(f)

    except Exception:
        return {}

    conditions: dict[int, object] = {}

    for comment in cm_tree.findall(f".//{{{_NS}}}comment"):
        ref = comment.get("ref", "")
        parts = [t.text for t in comment.findall(f".//{{{_NS}}}t") if t.text]
        text  = "".join(parts)

        m = re.search(
            r"Requirement level\s*:\s*Mandatory\s*\((.+?)\)",
            text, re.IGNORECASE | re.DOTALL
        )
        if not m:
            continue

        cond_text = m.group(1).strip()
        cidx      = _col_idx(ref)

        # "If samp_category = X"  — mandatory only when samp_category matches
        if_m = re.match(r"if\s+samp_category\s*=\s*(.+)", cond_text, re.IGNORECASE)
        if if_m:
            required = if_m.group(1).strip().rstrip(")").strip().lower()
            conditions[cidx] = lambda cat, rc=required: cat.strip().lower() == rc
            continue

        # "Mandatory unless samp_category = X, Y, or Z"
        unless_m = re.search(r"unless\s+samp_category\s*=\s*(.+)", cond_text, re.IGNORECASE)
        if unless_m:
            cats_str = unless_m.group(1)
            exempt   = {
                c.strip().strip(")").lower()
                for c in re.split(r",\s*(?:or\s+)?|\bor\b", cats_str)
                if c.strip()
            }
            conditions[cidx] = lambda cat, ec=exempt: cat.strip().lower() not in ec

    return conditions


# Tabs where rows represent samples/taxa and mandatory-column checks make sense.
# README and projectMetadata have non-standard layouts; otuRaw/otuFinal have
# samples as columns — these are handled separately.
_FAIRE_ROW_TABS = {"sampleMetadata", "experimentRunMetadata", "taxaRaw", "taxaFinal"}


def _check_mandatory_columns(rows: list, tab: str, fname: str,
                              results: Results, assay: str, d: str,
                              conditions: dict | None = None):
    """
    Parse the FAIRe requirements row (first cell = '# requirement_level_code',
    remaining cells = M / HR / R / O) to find mandatory columns, then verify
    every data row has a value in those columns.

    `conditions` maps col_index -> condition_fn(samp_category: str) -> bool.
    When a condition is present, the column is only mandatory for rows where
    condition_fn returns True (e.g. neg_cont_type is only M for negative controls).

    Control/negative samples (matched by CONTROL_PATTERNS or samp_category) may
    have blank mandatory fields — these are noted as PASS rather than flagged.

    FAIRe row layout:
        Row with first cell '# requirement_level_code': requirement codes per column
        Other '# ...' rows: section labels, descriptions, etc.
        First non-'#' row: column names (header)
        Remaining rows: data
    """
    _REQ_LEVELS = {"M", "HR", "R", "O"}
    if conditions is None:
        conditions = {}

    req_row        = None
    header_row     = None
    header_row_idx = None

    for i, row in enumerate(rows):
        first = str(row[0]).strip() if row[0] is not None else ""
        if first.startswith("#"):
            non_first = [str(c).strip() for c in row[1:] if c is not None and str(c).strip()]
            if non_first and all(v in _REQ_LEVELS for v in non_first):
                req_row = row
            continue
        if first == "":
            continue
        header_row     = row
        header_row_idx = i
        break

    if req_row is None or header_row is None or header_row_idx is None:
        results.add(assay, d, WARN,
            f"[{fname}] '{tab}': could not locate requirements row — "
            f"skipping mandatory-column check")
        return

    # Collect indices of mandatory (M) columns
    mandatory_indices = [
        idx for idx, v in enumerate(req_row)
        if v is not None and str(v).strip() == "M"
    ]
    if not mandatory_indices:
        return

    mandatory_names = [
        str(header_row[idx]).strip()
        if idx < len(header_row) and header_row[idx] is not None
        else f"col{idx}"
        for idx in mandatory_indices
    ]

    # Locate samp_name and samp_category columns
    header_lower = [str(v).lower().strip() if v is not None else "" for v in header_row]
    samp_name_idx     = header_lower.index("samp_name")     if "samp_name"     in header_lower else None
    samp_cat_idx      = header_lower.index("samp_category") if "samp_category" in header_lower else None

    _CONTROL_CATS = {"negative control", "positive control", "pcr standard"}

    # Scan every data row
    missing_by_col: dict[str, list[str]] = {}   # col_name -> [row_ids with missing value]
    control_blanks: dict[str, list[str]] = {}   # sample_name -> [col_names that are blank]

    for row in rows[header_row_idx + 1:]:
        if not row or all(v is None for v in row):
            continue

        # Determine sample name and category for this row
        samp_name = ""
        samp_cat  = ""
        if samp_name_idx is not None and samp_name_idx < len(row):
            samp_name = str(row[samp_name_idx]).strip() if row[samp_name_idx] is not None else ""
        if samp_cat_idx is not None and samp_cat_idx < len(row):
            samp_cat  = str(row[samp_cat_idx]).strip()  if row[samp_cat_idx]  is not None else ""

        # A row is a control if samp_category says so OR if the name contains a control pattern
        is_control = (
            samp_cat.lower() in _CONTROL_CATS
            or any(pat in samp_name for pat in CONTROL_PATTERNS)
        )

        row_id = samp_name or (str(row[0]).strip() if row[0] is not None else "(blank)")

        for col_idx, col_name in zip(mandatory_indices, mandatory_names):
            # If this column has a conditional rule, check whether it applies to this row
            cond_fn = conditions.get(col_idx)
            if cond_fn is not None and not cond_fn(samp_cat):
                continue  # Column not mandatory for this row's samp_category

            val   = row[col_idx] if col_idx < len(row) else None
            blank = (val is None
                     or str(val).strip() == ""
                     or str(val).strip().lower() in ("na", "n/a"))
            if blank:
                if is_control:
                    control_blanks.setdefault(row_id, []).append(col_name)
                else:
                    if col_name != "pcr_plate_id":
                        missing_by_col.setdefault(col_name, []).append(row_id)

    # Report failures
    if missing_by_col:
        for col_name, row_ids in list(missing_by_col.items()):
            results.add(assay, d, FAIL,
                f"[{fname}] '{tab}': mandatory column '{col_name}' missing values "
                f"in {len(row_ids)} row(s): {row_ids}")
    else:
        results.add(assay, d, PASS,
            f"[{fname}] '{tab}': all mandatory columns populated")

    # Note control-sample blanks as expected behaviour
    if control_blanks:
        n           = len(control_blanks)
        sample_list = sorted(control_blanks.keys())
        results.add(assay, d, PASS,
            f"[{fname}] '{tab}': {n} control/negative sample(s) have blank mandatory "
            f"fields — expected (controls may lack sample-collection metadata): "
            f"{sample_list}")


# ── 07-multiqc ─────────────────────────────────────────────────────────────────

def check_multiqc_dir(mq_dir: Path, p: str, a: str, results: Results):
    d = "07-multiqc"
    check_file(mq_dir / f"{p}_{a}_multiqc_report.html", results, a, d)
    #check_dir(mq_dir / f"{p}_{a}_multiqc_data",  results, a, d)
    #check_dir(mq_dir / f"{p}_{a}_multiqc_plots", results, a, d)


# ── 07-proportional_filter ─────────────────────────────────────────────────────

def check_proportional_filter_dir(pf_dir: Path, p: str, a: str, results: Results):
    d = "07-proportional_filter"
    exp = shared_expected(p, a)["07-proportional_filter"]

    for fname in exp:
        fpath = pf_dir / fname
        if not check_file(fpath, results, a, d):
            continue

        try:
            if "_proportional_stats.txt" in fname:
                text = fpath.read_text()
                if "Number of ASVs" not in text:
                    results.add(a, d, WARN, f"Unexpected content format in {fname}")
                else:
                    results.add(a, d, PASS, f"Stats content OK in {fname}")

            elif "_OTU_filtered.tsv" in fname:
                df = pd.read_csv(fpath, sep="\t", nrows=1)
                dupes = df.columns[df.columns.duplicated()].tolist()
                if dupes:
                    results.add(a, d, FAIL,
                        f"Duplicate sample columns in {fname}: {dupes}")
                else:
                    results.add(a, d, PASS, f"Sample columns unique in {fname}")
        except Exception as e:
            results.add(a, d, WARN, f"Could not parse {fname}: {e}")

    flag_unexpected(pf_dir, set(exp), results, a, d)


# ── --all only checks ──────────────────────────────────────────────────────────

def check_cutadapt_dir(ca_dir: Path, p: str, a: str, results: Results):
    d = "01-cutadapt"
    for subdir in ["assigned", "all-primers-trimmed", "unknown"]:
        check_dir(ca_dir / subdir, results, a, d)


def check_fastqc_dir(fq_dir: Path, p: str, a: str, results: Results):
    d = "01-fastqc"
    check_dir(fq_dir, results, a, d)
    if fq_dir.exists():
        htmls = list(fq_dir.glob("*.html"))
        zips  = list(fq_dir.glob("*.zip"))
        if not htmls:
            results.add(a, d, FAIL, "No FastQC .html files found")
        else:
            results.add(a, d, PASS, f"{len(htmls)} FastQC .html files found")
        if not zips:
            results.add(a, d, WARN, "No FastQC .zip files found")
        else:
            results.add(a, d, PASS, f"{len(zips)} FastQC .zip files found")


def check_seqkit_dir(sq_dir: Path, p: str, a: str, results: Results):
    d = "01-seqkit_stats"
    exp = all_only_expected(p, a)["01-seqkit_stats"]
    for fname in exp:
        check_file(sq_dir / fname, results, a, d)
    flag_unexpected(sq_dir, set(exp), results, a, d)


def check_dada2_dir(dada_dir: Path, p: str, a: str, results: Results):
    d = "02-dada2"
    for fname in all_only_expected(p, a)["02-dada2"]:
        check_file(dada_dir / fname, results, a, d)
    check_dir(dada_dir / "plots", results, a, d)


def check_lulu_dir(lulu_dir: Path, p: str, a: str, results: Results):
    d = "03-lulu"
    for fname in all_only_expected(p, a)["03-lulu"]:
        check_file(lulu_dir / fname, results, a, d)
    check_dir(lulu_dir / f"{p}_{a}_asv_db", results, a, d)


def check_blast_dir(blast_dir: Path, p: str, a: str, results: Results):
    d = "04-blast"
    for fname in all_only_expected(p, a)["04-blast"]:
        check_file(blast_dir / fname, results, a, d)


def check_ocomnbc_dir(nbc_dir: Path, p: str, a: str, results: Results):
    d = "04-ocomnbc"
    for fname in all_only_expected(p, a)["04-ocomnbc"]:
        check_file(nbc_dir / fname, results, a, d)


def check_pipeline_info_dir(pi_dir: Path, p: str, a: str, results: Results):
    d = "07-pipeline_info"
    for fname in [f"{p}_{a}_samplesheet.valid.csv", "software_versions.yml"]:
        check_file(pi_dir / fname, results, a, d)
    for pattern, label in [
        ("execution_report_*.html",    "execution_report"),
        ("execution_timeline_*.html",  "execution_timeline"),
        ("execution_trace_*.txt",      "execution_trace"),
        ("pipeline_dag_*.html",        "pipeline_dag"),
    ]:
        matches = list(pi_dir.glob(pattern)) if pi_dir.exists() else []
        if not matches:
            results.add(a, d, FAIL, f"MISSING: no {label} file found")
        else:
            results.add(a, d, PASS, f"OK: {matches[0].name}")


# ── assay validator ────────────────────────────────────────────────────────────

def validate_assay(assay_dir: Path, project_id: str, check_all: bool,
                   results: Results, sing2: str, log: logging.Logger):
    a = assay_dir.name
    log.info(f"{'='*60}")
    log.info(f"Validating assay: {a}")
    log.info(f"{'='*60}")

    p = project_id

    # Always checked (shared with collaborators)
    check_lca_dir(assay_dir / "05-lca", p, a, results)
    check_aquamap_dir(assay_dir / "06-aquamap", p, a, results)
    check_phyloseq_dir(assay_dir / "06-phyloseq", p, a,
                       assay_dir / "07-faire", results, sing2)
    check_faire_dir(assay_dir / "07-faire", p, a, results)
    check_multiqc_dir(assay_dir / "07-multiqc", p, a, results)
    check_proportional_filter_dir(assay_dir / "07-proportional_filter", p, a, results)

    if check_all:
        check_cutadapt_dir(assay_dir / "01-cutadapt", p, a, results)
        check_fastqc_dir(assay_dir / "01-fastqc", p, a, results)
        check_seqkit_dir(assay_dir / "01-seqkit_stats", p, a, results)
        check_dada2_dir(assay_dir / "02-dada2", p, a, results)
        check_lulu_dir(assay_dir / "03-lulu", p, a, results)
        check_blast_dir(assay_dir / "04-blast", p, a, results)
        check_ocomnbc_dir(assay_dir / "04-ocomnbc", p, a, results)
        check_pipeline_info_dir(assay_dir / "07-pipeline_info", p, a, results)

    # Per-assay mini-summary
    items  = results.for_assay(a)
    n_fail = sum(1 for _, _, s, _ in items if s == FAIL)
    n_warn = sum(1 for _, _, s, _ in items if s == WARN)
    n_pass = sum(1 for _, _, s, _ in items if s == PASS)
    log.info(f"Assay {a} — {n_pass} PASS  {n_warn} WARN  {n_fail} FAIL")


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="OceanOmics eDNA Pipeline Output Validator (OcOm_B218)"
    )
    parser.add_argument(
        "project_dir", type=Path,
        help="Project directory containing assay subdirectories (e.g. /path/to/OcOm_2516)"
    )
    parser.add_argument(
        "--partial", dest="partial", action="store_true",
        help="Check shared/collaborator dirs only (default: Check all dirs)"
    )
    parser.add_argument(
        "--log", type=Path,
        help="Log file path (default: <project_dir>/validation_<timestamp>.log)"
    )
    parser.add_argument(
        "--sing2", type=str,
        default=os.environ.get("SING2", ""),
        help="Path to singularity images directory (default: $SING2)"
    )
    args = parser.parse_args()

    if args.sing2 == "":
        args.sing2 = os.environ.get("SING", "")

    project_dir = args.project_dir.resolve()
    if not project_dir.is_dir():
        print(f"ERROR: {project_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    project_id = str(project_dir).split("/")[-1]
    timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path   = args.log or project_dir / f"validation_{project_id}_{timestamp}.log"

    log = setup_logging(log_path)
    log.info("OceanOmics eDNA Output Validator  |  SOP OcOm_B218")
    log.info(f"Project   : {project_id}")
    log.info(f"Mode      : {'shared directories only' if args.partial else 'all directories'}")
    log.info(f"SING2     : {args.sing2 or '(not set — phyloseq .rds checks will be skipped)'}")
    log.info(f"Log file  : {log_path}")

    # Discover assay directories
    assay_dirs = sorted(
        d for d in project_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    )
    assay_dirs = []
    for d in project_dir.iterdir():
        path = Path(d)
        #for entry in path.rglob("*"):
        #print(path.iterdir())
        subdirs = [x for x in Path(d).rglob('*') if x.is_dir()]
        for i in subdirs:
            assay_dirs.append(i)
    if not assay_dirs:
        log.error(f"No subdirectories found in {project_dir}")
        sys.exit(1)
    
    # Filter to dirs that look like assay outputs
    assay_dirs = [
        d for d in assay_dirs
        if any((d / sub).exists() for sub in ("05-lca", "06-phyloseq", "07-faire"))
    ]

    if not assay_dirs:
        log.error("No assay directories found (expected subdirs: 05-lca, 06-phyloseq, 07-faire)")
        sys.exit(1)

    log.info(f"Assays    : {[d.name for d in assay_dirs]}")
    log.info(f"Run date  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    directories = [x.name for x in project_dir.iterdir() if x.is_dir()]
    if "logs" in directories:
        directories.remove("logs")
    log.info(f"Directory : {str(directories)}")

    results = Results(log)

    check_all = True
    if args.partial:
        check_all = False

    for assay_dir in assay_dirs:
        validate_assay(assay_dir, project_id, check_all, results, args.sing2, log)

    # ── final summary ──────────────────────────────────────────────────────────
    fails, warns, passes = results.counts()
    log.info("")
    log.info("=" * 60)
    log.info("FINAL SUMMARY")
    log.info("=" * 60)
    log.info(f"Total checks : {len(results.all_items())}")
    log.info(f"PASS         : {len(passes)}")
    log.info(f"WARN         : {len(warns)}")
    log.info(f"FAIL         : {len(fails)}")
    log.info(f"Run date     : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log.info(f"Directory    : {str(directories)}")

    if fails:
        log.error("")
        log.error(f"FAILURES ({len(fails)}):")
        for assay, directory, msg in fails:
            log.error(f"  [{assay}] {directory}: {msg}")

    if warns:
        log.warning("")
        log.warning(f"WARNINGS ({len(warns)}):")
        for assay, directory, msg in warns:
            log.warning(f"  [{assay}] {directory}: {msg}")

    log.info("")
    if fails:
        log.error("RESULT: FAIL — data is NOT ready to share")
    elif warns:
        log.info("RESULT: PASS WITH WARNINGS — review warnings before sharing")
    else:
        log.info("RESULT: ALL CHECKS PASSED — data is ready to share")

    log.info(f"\nLog written to: {log_path}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
