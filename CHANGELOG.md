# mskcc/neoantigenpipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.4.0 - [date]

### `Added`

- New `GENERATE_MUTATED_PEPTIDES` subworkflow that groups `MUTALYZER_RETRIEVER`, `GENERATEMUTFASTA`, `GENERATEHLASTRING`, and `NEOSV` into a single reusable unit, replacing individual module calls in the main workflow ([#82](https://github.com/mskcc/neoantigen-pipeline/pull/82), [#85](https://github.com/mskcc/neoantigen-pipeline/pull/85))
- New `MUTALYZER_RETRIEVER` module for retrieving and caching Mutalyzer annotation results; controlled by new `mutalyzer_cache_dir` parameter ([#85](https://github.com/mskcc/neoantigen-pipeline/pull/85))
- New `NEOSV` module for structural variant neoantigen FASTA generation ([#85](https://github.com/mskcc/neoantigen-pipeline/pull/85))
- New `run_phylowgs` parameter (boolean) to make PhyloWGS tumor phylogeny inference optional; when disabled, a lightweight stub replaces PhyloWGS ([#80](https://github.com/mskcc/neoantigen-pipeline/pull/80))
- New `PHYLOWGS_STUB` local module to produce placeholder PhyloWGS outputs when `run_phylowgs = false` ([#80](https://github.com/mskcc/neoantigen-pipeline/pull/80))
- New `conf/juno.config` production profile for the MSKCC Juno HPC cluster (LSF executor, Singularity, cluster-specific reference paths)
- New parameters `reference_fasta` and `reference_gff3` for Mutalyzer genome annotation input
- New `pipeline_options` section in `nextflow_schema.json` documenting all pipeline-specific parameters with descriptions
- End-to-end pipeline test (`tests/default.nf.test`) with snapshot validation
- `.nftignore` rules for non-deterministic pipeline outputs
- `CLAUDE.md` developer documentation with architecture overview, commands, and key parameter reference
- nf-core template synced from v3.2.0 to v3.5.2 — includes new GitHub Actions for sharded nf-test runs, `nf-test.yml` workflow, and `template-version-comment.yml`

### `Changed`

- `GENERATE_MUTATED_PEPTIDES` subworkflow now orchestrates mutated FASTA generation; the `NETMHCSTABANDPAN` subworkflow interface updated accordingly to accept pre-generated FASTA + HLA channels
- `generatemutfasta` module relocated from `modules/msk/neoantigenutils/generatemutfasta/` to `modules/msk/generatemutfasta/` (directory flattened to satisfy nf-core linting)
- `prod.config` updated: replaced `cds` reference with `reference_fasta` + `reference_gff3`; `run_phylowgs` defaults to `false`; `config_profile_description` corrected to "Production conf to run the pipeline"
- Workflow-level `new Tuple(...)` calls replaced with list literals `[...]` for Nextflow ≥ 25.04.0 compatibility
- Local variables in `merge_for_input_generation()` now use explicit `def` keywords as required by Nextflow ≥ 25.04.0
- Software version collection migrated to use `Channel.topic("versions")` for broader process coverage
- `nf-core lint` CI switched from SSH to HTTPS for modules repository URL
- Bumped minimum Nextflow version requirement

### `Fixed`

- Heredoc `END_VERSIONS` blocks: added tab indentation and removed empty lines to resolve Nextflow ≥ 25.04.0 heredoc parsing errors
- Trailing whitespace stripped from heredoc padding lines
- `nf-test` snapshot mismatches resolved after module and Nextflow version upgrades
- README formatting fixes and added tip for running pipeline with the IRIS Nextflow config
- Duplicate parameter definitions removed from config files
- `fix-linting.yml` workflow renamed to `fix_linting.yml` for consistency

### `Dependencies`

| Module / Tool       | Change                          |
| ------------------- | ------------------------------- |
| `mutalyzer`         | Added (new module)              |
| `neosv`             | Added (new module)              |
| `generatemutfasta`  | Relocated; script updated       |
| `multiqc`           | Updated to nf-core latest       |
| `nf-schema`         | v2.3.0 (via nf-core template)   |
| Nextflow            | Minimum version bumped          |
| nf-core template    | Synced 3.2.0 → 3.5.2            |

### `Deprecated`

- `cds` parameter removed from production configs; replaced by `reference_fasta` + `reference_gff3`
- `modules/msk/neoantigenutils/generatemutfasta/` path deprecated; use `modules/msk/generatemutfasta/`

---

## v1.3.1 - 2025-03-27

### `Added`

- Pipeline-level parameters `compute_fitness_a`, `compute_fitness_k`, `compute_fitness_w`, `kd_cutoff`, `gtf`, `cdna`, and `cds` added to `prod.config` and `nextflow_schema.json` for explicit schema validation and documentation
- `trace_report_suffix` schema parameter added (nf-core template)

### `Fixed`

- `neoantigeninput` updated to correctly handle edge cases in neoantigen input generation
- Corrected the `generate_input.py` script when processing multi-sample inputs (`updated generateinput`)

### `Dependencies`

| Module / Tool          | Change                        |
| ---------------------- | ----------------------------- |
| `neoantigeninput`      | Script updated (v1.9)         |

---

## v1.3.0 - 2025-03-10

### `Added`

- `gtf` and `cdna` reference parameters passed to `NEOANTIGENUTILS_NEOANTIGENINPUT` to enable transcript-level annotation during neoantigen input generation
- `compute_fitness_a`, `compute_fitness_k`, `compute_fitness_w` parameters added to `prod.config`

### `Changed`

- Replaced standalone `NETMHCPAN` module call with the combined `NETMHCSTABANDPAN` subworkflow for all binding prediction steps
- `neoantigeninput` updated to use `pyensembl`-backed gene/transcript annotation (`generate_input.py` v1.9)

### `Deprecated`

- `NETMHCPAN` standalone module removed from workflow

### `Dependencies`

| Module / Tool          | Change                                       |
| ---------------------- | -------------------------------------------- |
| `netmhcpan`            | Removed as standalone module                 |
| `neoantigeninput`      | Updated with pyensembl transcript annotation |
| `computefitness`       | Updated                                      |
| `convertannotjson`     | Updated                                      |
| `formatnetmhcpan`      | Updated                                      |
| `generatehlastring`    | Updated                                      |
| `generatemutfasta`     | Updated                                      |

---

## v1.2.0 - 2025-02-13

### `Added`

- `kd_cutoff` parameter added to control binding affinity filtering threshold; wired into `NEOANTIGENUTILS_NEOANTIGENINPUT` via `ext.args` in `modules.config`

### `Fixed`

- `generate_input.py` (v1.9): refactored internal variable naming (`noposID` → `no_positon_ID`) and fixed handling of no-position identity lookups in WT peptide dictionaries

### `Dependencies`

| Module / Tool          | Change                                      |
| ---------------------- | ------------------------------------------- |
| `neoantigeninput`      | Script fix (v1.8 → v1.9)                   |
| `multiqc`              | Updated to nf-core latest                   |
| nf-core template       | Synced to 3.2.0                             |

---

## v1.1.1 - 2025-01-29

### `Added`

- `phylo_num_chains` parameter added to `prod.config`
- `netmhc3` flag added to `prod.config` to enable NetMHC 3-column output format by default in production

### `Fixed`

- Versions output file renamed from generic `pipeline_software_mqc_versions.yml` to `neoantigenpipeline_software_mqc_versions.yml`
- Singularity profile: disabled Apptainer fallback to fix container engine conflicts in CI
- Reverted unstable `download_pipeline` config change that broke CI

### `Dependencies`

| Module / Tool     | Change                                       |
| ----------------- | -------------------------------------------- |
| nf-core template  | Synced through 3.1.2                         |

---

## v1.1.0 - 2024-12-06

### `Added`

- `NEOANTIGENUTILS_CONVERTANNOTJSON` module added to produce a user-friendly TSV summary alongside the annotated JSON output
- `netmhc3` parameter added to support NetMHC 3.4 output format parsing in `formatnetmhcpan`
- NetMHC 3.4 (`netmhc3`) handling added to the `NETMHCSTABANDPAN` subworkflow
- Workflow diagram updated in README
- `phylo_num_chains` set to 15 in `prod.config`

### `Changed`

- `prod.config` migrated from legacy `max_cpus` / `max_memory` / `max_time` params to `resourceLimits` block (nf-core template 3.x)
- JSON schema upgraded from draft-07 to draft 2020-12; `definitions` key replaced with `$defs`
- Removed `max_job_request_options` schema section (superseded by `resourceLimits`)
- `nf-validation` plugin replaced with `nf-schema`
- Switched output channel metadata key from `typePan` to `fromPan`

### `Fixed`

- Empty SV channel handling: properly initialise empty BEDPE input for `NEOANTIGENUTILS_NEOANTIGENINPUT`
- `neoantigeninput` updated to handle optional BEDPE input and edge cases
- Various nf-core linting and Prettier formatting fixes

### `Dependencies`

| Module / Tool          | Change                             |
| ---------------------- | ---------------------------------- |
| `convertannotjson`     | Added (new module)                 |
| `nf-schema`            | Replaces `nf-validation`           |
| nf-core template       | Synced through 3.1.x               |

---

## v1.0.0 - 2024-07-31

Initial release of mskcc/neoantigenpipeline, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- End-to-end neoantigen quality pipeline implementing the Luksza et al. neoantigen editing and fitness framework
- **PhyloWGS** subworkflow: `PARSECNVS` → `CREATEINPUT` → `MULTIEVOLVE` → `WRITERESULTS` for tumor subclonal reconstruction from MAF + Facets CNCF input
- **NetMHCstabandpan** subworkflow: generates mutant FASTA from MAF + CDS/cDNA references, runs `netMHCpan-4.1` (binding affinity) and `netMHCstabpan` (stability), producing MUT and WT predictions
- **NeoantigenInput** step: merges PhyloWGS clonal tree with netMHCpan results via `NEOANTIGENUTILS_NEOANTIGENINPUT`
- **NeoantigenEditing** subworkflow: `ALIGNTOIEDB` → `COMPUTEFITNESS` → `CONVERTANNOTJSON`, computing R, logC, logA, and quality scores
- CSV samplesheet input format: `sample`, `maf`, `facets_hisens_cncf`, `hla_file`
- `prod.config` for MSKCC production (LSF executor, Singularity, GRCh37 reference paths)
- `test.config` for CI with reduced PhyloWGS iterations
- nf-schema parameter validation (`nextflow_schema.json`)
- MultiQC and software version reporting
