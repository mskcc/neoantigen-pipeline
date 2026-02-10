# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

mskcc/neoantigenpipeline is a **Nextflow DSL2** bioinformatics pipeline (v1.3.1) that computes neoantigen qualities from DNA and RNA-Seq data. Built on the nf-core template (v3.2.0), it implements Luksza et al.'s neoantigen editing and fitness framework for MSKCC investigators. Input is a CSV samplesheet with MAF, Facets CNCF, and HLA files; output is annotated JSON with clonal structure, neoantigen burden, and fitness scores plus a TSV summary.

## Commands

```bash
# Run pipeline (test)
nextflow run . -profile test,docker --outdir ./results

# Run pipeline (production - LSF executor at MSKCC)
nextflow run . -profile prod,singularity --input samplesheet.csv --outdir ./results

# Run all nf-tests
nf-test run

# Run a single module test
nf-test run modules/msk/phylowgs/createinput/tests/main.nf.test

# Run a subworkflow test
nf-test run subworkflows/msk/phylowgs/tests/main.nf.test

# Lint (nf-core)
nf-core pipelines lint

# Pre-commit hooks (Prettier + EditorConfig)
pre-commit run --all-files
```

## Architecture

### Pipeline Stages (in order)

The main workflow is in `workflows/neoantigenpipeline.nf`. Data flows through four major stages:

1. **PhyloWGS** — Reconstructs subclonal composition from MAF + Facets copy number data
   - `parsecnvs` → `createinput` → `multievolve` → `writeresults`
   - Outputs: summary JSON, mutations JSON, mutation assignments

2. **netMHCstabandpan** — Predicts MHC peptide binding (subworkflow combines two tools)
   - Generates mutant FASTA from MAF + CDS/cDNA references
   - Runs netMHCpan-4.1 for binding affinity (Kd) and netMHCstabpan for stability
   - Produces MUT and WT predictions, filtered by `meta.typeMut` and `meta.fromStab` flags

3. **NeoantigenInput** — Merges PhyloWGS tree with netMHCpan results
   - `neoantigenutils/neoantigeninput` joins all channels by sample ID using `merge_for_input_generation()`
   - Applies `kd_cutoff` parameter to filter weak binders

4. **NeoantigenEditing** — Computes fitness and quality scores
   - `aligntoiedb` aligns peptides against IEDB database
   - `computefitness` calculates R, logC, logA, and quality scores (params: `a`, `k`, `w`)
   - `convertannotjson` produces final TSV output

### Module Organization

```
modules/msk/{tool}/{process}/main.nf    # Custom MSK modules
modules/nf-core/multiqc/                 # nf-core standard modules
subworkflows/msk/{name}/main.nf          # MSK subworkflows combining modules
subworkflows/nf-core/                    # nf-core utility subworkflows
```

Each module has a `tests/` directory with `main.nf.test` files for nf-test.

### Channel Merging Pattern

The workflow uses a custom `merge_for_input_generation()` function that joins six channels by `meta.id` to synchronize PhyloWGS outputs with netMHCpan results before neoantigen input generation. This is the most complex data-flow logic in the pipeline.

### Key Configuration

- `conf/base.config` — Default resource labels (`process_single`, `process_low`, `process_medium`, etc.)
- `conf/prod.config` — MSKCC production: LSF executor, singularity, reference genome URLs, PhyloWGS MCMC params, fitness model params
- `conf/test.config` — Minimal test: reduced PhyloWGS iterations (2 burnin, 2 MCMC, 2 chains)
- `conf/modules.config` — Per-process `ext.args` (e.g., facets format for parsecnvs, MCMC params for multievolve)
- `nextflow_schema.json` — Parameter validation schema (nf-schema plugin v2.3.0)

### Key Parameters

| Parameter                                                        | Purpose                                      |
| ---------------------------------------------------------------- | -------------------------------------------- |
| `phylo_burnin_samples`, `phylo_mcmc_samples`, `phylo_num_chains` | PhyloWGS MCMC settings                       |
| `kd_cutoff`                                                      | Binding affinity threshold (default: 500 nM) |
| `compute_fitness_a`, `compute_fitness_k`, `compute_fitness_w`    | Fitness model parameters                     |
| `iedbfasta`, `cds`, `cdna`, `gtf`                                | Reference data URLs                          |

### Input Samplesheet Format

```csv
sample,maf,facets_hisens_cncf,hla_file
tumor_normal,mutations.maf,facets_hisens.cncf.txt,winners.hla.txt
```

## CI/CD

- **CI** (`.github/workflows/ci.yml`): Runs pipeline with test profile on push/PR to dev/main/master. Tests Nextflow 24.04.0 and latest.
- **Linting** (`.github/workflows/linting.yml`): Pre-commit hooks + `nf-core pipelines lint`. Uses `--release` mode on master.
- **Pre-commit**: Prettier v3.1.0 and EditorConfig checker v3.1.2.

## Nextflow Version

Requires Nextflow >= 24.04.0. Uses DSL2 with `nf-schema@2.3.0` plugin.
