# SOP OcOm_B218 — eDNA Pipeline Output Validation

**Purpose:** Audit eDNA metabarcoding pipeline outputs before sharing data with collaborators  
**Script:** `validate_edna_output.py`  
**Helper:** `check_phyloseq.R`  
**Platform:** Setonix HPC, Pawsey Supercomputing Centre

---

## 1. Prerequisites

### 1.1 Software

`singularity` is loaded automatically at Setonix login. `uv` must be installed once per user — it is not a Setonix module.

**Install `uv` (one-time, run on a login node):**

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

The installer places `uv` at `~/.local/bin/uv` and adds `~/.local/bin` to your `PATH` in `~/.bashrc`. Apply immediately with `source ~/.bashrc`.

| Tool | Notes |
|------|-------|
| `uv` | Inline script runner — handles Python dependencies automatically. Install as above. |
| `singularity` | Required for phyloseq `.rds` checks — loaded at Setonix login |
| Python ≥ 3.11 | Managed by `uv` — no separate install needed |

### 1.2 Script files

Both files must be in the same directory:

```
validate_edna_output.py
check_phyloseq.R
```

### 1.3 Environment variables

Add the following to your `~/.bashrc` on Setonix. Both variables must be set for all checks to run correctly.

```bash
# Redirect uv's package cache to scratch — home directory quota is too small
export UV_CACHE_DIR=/scratch/pawsey0964/$USER/.cache/uv

# Path to shared Singularity images — required for phyloseq .rds checks
# This is already set at Setonix login; add the line below only if it is missing
export SING=/software/projects/pawsey0964/singularity
```

Apply changes immediately:

```bash
source ~/.bashrc
```

> **Note:** If `UV_CACHE_DIR` is not set, `uv` will write to `~/.cache/uv` and the run will fail when the home directory quota is exceeded.  
> **Note:** If `SING` is not set, phyloseq `.rds` checks are skipped with a WARN.

---

## 2. Running the Validator

### 2.1 Basic usage

```bash
# Check all pipeline directories including intermediate steps
uv run --script validate_edna_output.py <project_dir>

# Check shared/collaborator directories only (recommended before sharing)
uv run --script validate_edna_output.py <project_dir> --partial

# Next, run backup script
bash backup_results.sh <project_dir>
```

**Example:**

```bash
uv run --script validate_edna_output.py <project_dir>
```

### 2.2 Optional arguments

```bash
# Write log to a custom path instead of the default location
uv run --script validate_edna_output.py <project_dir> --log /path/to/output.log

# Override $SING to specify a different Singularity images directory
uv run --script validate_edna_output.py <project_dir> --sing /path/to/sifs
```

### 2.3 Output

The log is written to:

```
<project_dir>/validation_<project_id>_<timestamp>.log
```

The final line of the log gives an overall result:

| Result | Meaning |
|--------|---------|
| `RESULT: ALL CHECKS PASSED — data is ready to share` | No issues found |
| `RESULT: PASS WITH WARNINGS — review warnings before sharing` | Minor issues; review before sharing |
| `RESULT: FAIL — data is NOT ready to share` | Critical issues must be resolved first |

---

## 3. Check Modes

| Mode | What is checked |
|------|----------------|
| `--partial` | Shared/collaborator directories: `05-lca`, `06-aquamap`, `06-phyloseq`, `07-faire`, `07-multiqc`, `07-proportional_filter` |
| Default | All of the above **plus** early pipeline directories: `01-cutadapt`, `01-fastqc`, `01-seqkit_stats`, `02-dada2`, `03-lulu`, `04-blast`, `04-ocomnbc`, `07-pipeline_info` |

---

## 4. Result Levels

| Level | Meaning |
|-------|---------|
| `PASS` | Check passed |
| `WARN` | Potential issue; review before sharing but does not block |
| `FAIL` | Critical issue; data must NOT be shared until resolved |

---

## 5. Checks Performed

Assay directories are auto-discovered: any subdirectory of the project root containing `05-lca`, `06-phyloseq`, or `07-faire` is treated as an assay (e.g. `16SFish`, `MiFishU`, `COILeray`).

---

### Partial Mode — Shared / Collaborator Directories

#### `05-lca` — LCA taxonomy files

- All expected TSV files are present and non-empty
- `_taxa_final.tsv` files: `seq_id` column is present and values are unique (no duplicate ASVs)
- `_lca_with_fishbase_output.tsv` files: all required taxonomy columns are present (`domain`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `otu`, `querycoverage`, `%id`)
- Unexpected files flagged as WARN

#### `06-aquamap` — AquaMaps species range CSVs

- All expected CSV files are present and non-empty
- Species index (row names) are unique — no duplicate species entries
- All-NA values flagged as WARN (possible missing lat/lon data)
- Files containing `no_aquamap` in the name flagged as WARN (missing coordinate data)
- Unexpected files flagged as WARN

#### `06-phyloseq` — phyloseq R objects and taxonomy TSVs

- All expected `.rds` and `_final_taxa.tsv` files are present and non-empty
- `_final_taxa.tsv` files: required taxonomy columns present (`domain`, `phylum`, `class`, `order`, `family`, `genus`, `species`, `lca`)
- Each `_flagged_phyloseq.rds` file is loaded and checked via R (Singularity container):
  - `otu_table`, `sample_data`, and `tax_table` components are present
  - Sample names in `otu_table` are unique
- Sample name cross-check (phyloseq `sam_data` vs. FAIRe `sampleMetadata`):
  - Control samples (`_NTC_`, `_ITC_`, `_EB_`, `_WC_`, `_FC_`) in phyloseq but absent from FAIRe → PASS (controls are excluded from FAIRe by design)
  - Non-control samples in phyloseq but absent from FAIRe → WARN
  - Samples in FAIRe but absent from phyloseq are classified by cause:
    - Non-ITC control samples (EB, WC, FC, NTC): PASS — controls may have 0 reads
    - ITC positive controls: WARN — ITC should always have reads; absence is concerning
    - Samples not present in the pipeline samplesheet: PASS — excluded/discarded before pipeline ran
    - Samples in the pipeline samplesheet with `discarded=False`: WARN — undiscarded sample missing from phyloseq; check read counts
- Column name cross-check (`sam_data` vs. FAIRe `sampleMetadata`): mismatches flagged as WARN; internal pipeline columns (`discarded`, `sample_type`, `use_for_filter`) are excluded
- Unexpected files flagged as WARN

#### `07-faire` — FAIRe metadata xlsx files

- All expected `.xlsx` files are present and non-empty
- Required tabs present: `README`, `projectMetadata`, `sampleMetadata`, `experimentRunMetadata`, `taxaRaw`, `taxaFinal`, `otuRaw`, `otuFinal`
- Unexpected tabs flagged as WARN (known extras such as `Drop-down values` are not flagged)
- `projectMetadata`: project ID found in the values column
- Each tab is non-empty
- **Mandatory column completeness** — checked across `sampleMetadata`, `experimentRunMetadata`, `taxaRaw`, and `taxaFinal`:
  - The requirement level row (`# requirement_level_code`) is parsed to identify mandatory (`M`) columns
  - Conditional mandatory rules are read from the cell comments in the xlsx template, e.g.:
    - `neg_cont_type` is only mandatory when `samp_category = negative control`
    - `pos_cont_type` is only mandatory when `samp_category = positive control`
    - `decimalLatitude`, `decimalLongitude`, `env_broad_scale`, `env_local_scale`, `env_medium` are mandatory unless `samp_category` is a control type
  - Control/negative samples with blank mandatory fields → PASS with note (expected — controls lack sample-collection metadata)
  - Non-control rows with blank mandatory fields → FAIL
- `taxaFinal`: ASV IDs (first column) are unique
- `otuRaw` / `otuFinal`: sample column headers are unique (no duplicate sample names)
- Unexpected files flagged as WARN

#### `07-multiqc` — MultiQC report

- MultiQC HTML report file is present and non-empty
- `multiqc_data/` directory is present and non-empty
- `multiqc_plots/` directory is present and non-empty

#### `07-proportional_filter` — Proportional filter outputs

- All expected TSV/RDS/TXT files are present and non-empty
- `_proportional_stats.txt` files: content contains expected `Number of ASVs` summary text
- `_OTU_filtered.tsv` files: sample column headers are unique
- Unexpected files flagged as WARN

---

### Default Mode — Full Pipeline Directories

#### `01-cutadapt` — Primer trimming

- `assigned/`, `all-primers-trimmed/`, and `unknown/` subdirectories are present and non-empty

#### `01-fastqc` — Read quality reports

- Directory is present and non-empty
- At least one FastQC `.html` report is present
- FastQC `.zip` archives are present (WARN if missing)

#### `01-seqkit_stats` — Sequence statistics

- `assigned_seqkit_stats.txt`, `final_seqkit_stats.txt`, `raw_seqkit_stats.txt`, `unknown_seqkit_stats.txt`, and the project-specific prefilter stats file are present and non-empty
- Unexpected files flagged as WARN

#### `02-dada2` — DADA2 ASV calling

- Expected ASV FASTA, ASV table (TSV/CSV), LCA input TSV, and `.rds` sequence table are present and non-empty
- `<project>_<assay>_asv_track_reads.txt` is present and non-empty — checked at the top level of `02-dada2/` first, then inside `plots/` (location varies by pipeline version)
- `plots/` subdirectory is present and non-empty

#### `03-lulu` — LULU curation

- Expected curated ASV FASTA, curated table, LULU map, and match list files are present and non-empty
- BLAST database directory (`<project>_<assay>_asv_db/`) is present and non-empty

#### `04-blast` — BLAST results

- Raw and LULU-curated BLAST output files (both `blast` and `blast2` variants) are present and non-empty

#### `04-ocomnbc` — OceanOmics NBC classifier

- NBC output TSVs (all four blast/lulucurated variants) are present and non-empty

#### `07-pipeline_info` — Nextflow run info

- `<project>_<assay>_samplesheet.valid.csv` and `software_versions.yml` are present and non-empty
- Nextflow execution report, timeline, trace, and DAG files (wildcard-matched) are present
- The valid samplesheet is also used by the `06-phyloseq` check to classify samples that appear in FAIRe but are absent from the phyloseq object

---

## 6. Technical Notes

- The raw XML parser is used for all `.xlsx` reads to work around an openpyxl dimension bug present in FAIRe-format files. Cell comments in the xlsx are also parsed at runtime to extract conditional mandatory field rules — no hardcoding of conditions in the script.
- The valid samplesheet (`07-pipeline_info/<project>_<assay>_samplesheet.valid.csv`) contains only non-discarded samples that entered the pipeline. Sample names carry a `_T<N>` run-index suffix which is stripped before matching against FAIRe/phyloseq names.
