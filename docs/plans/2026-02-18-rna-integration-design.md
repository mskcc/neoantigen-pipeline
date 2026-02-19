# RNA-seq Integration Design for Neoantigen Pipeline

**Date**: 2026-02-18
**Branch**: feature/rna-integration-design
**Pipeline Version**: v1.4.0 (feature/module_updates)
**Status**: Design approved, pending implementation

## Overview

Integrate RNA-seq outputs from the FORTE pipeline into the neoantigen prediction pipeline to report RNA-sourced evidence alongside existing DNA-based neoantigen quality scores. The scoring model (Luksza et al.) remains unchanged; RNA data is annotation and reporting only.

## Requirements

- Do not modify the existing scoring model (R, Kd, logC, logA, quality)
- Report RNA data alongside final quality scores as additional columns
- Consume pre-computed FORTE output files (decoupled pipelines)
- Fully optional: pipeline works identically with DNA-only samples
- Include fusion neoantigens from FORTE as new predictions in a separate output file
- Produce both: key RNA columns in the main neoantigen TSV and a separate detailed RNA report

## Architecture

### Approach: Post-Processing Annotation (Approach A)

Selected over two alternatives:
- **Approach B** (inject RNA data at generate_input.py): Too invasive, modifies complex existing code with zero benefit since scoring is unchanged
- **Approach C** (hybrid early+late): Unnecessary complexity

The pipeline adds two new optional paths after the existing workflow:

```
Existing Pipeline (unchanged):
  MAF + HLA + Facets -> GENERATE_MUTATED_PEPTIDES -> NETMHCSTABANDPAN ->
  NEOANTIGENUTILS_NEOANTIGENINPUT -> NEOANTIGEN_EDITING -> CONVERTANNOTJSON
  -> neoantigen.tsv

New Path 1 - RNA Annotation (post-processing):
  neoantigen.tsv + MAF (with rna_* columns) + kallisto_abundance ->
  NEOANTIGEN_RNA_ANNOTATE -> annotated_neoantigen.tsv + rna_detail_report.tsv

New Path 2 - Fusion Neoantigens (parallel):
  agfusion_dir + HLA -> NEOANTIGEN_FUSION_PREPARE -> NETMHCSTABANDPAN ->
  NEOANTIGEN_EDITING -> fusion_neoantigen.tsv
```

## Samplesheet Changes

Three new optional columns added to `assets/schema_input.json`:

| Column | Required | Source | Format |
|--------|----------|--------|--------|
| `sample` | Yes | Existing | Sample ID |
| `maf` | Yes | Existing | MAF file (may contain FORTE rna_* columns) |
| `facets_hisens_cncf` | Yes | Existing | Facets CNV file |
| `hla_file` | Yes | Existing | Polysolver HLA output |
| `kallisto_abundance` | No (new) | FORTE Kallisto | `abundance.tsv` with transcript-level TPM |
| `agfusion_dir` | No (new) | FORTE AGFusion | Directory with fusion `.fa` files |

### RNA VAF: Auto-Detected from Input MAF

FORTE's GBCMS fillout workflow produces a combined MAF that merges RNA read counts into the original MAF. When the input MAF contains `rna_*` columns, the pipeline auto-detects and uses them. No separate fillout file is needed.

FORTE adds these columns via `rna_fillout_combine.R`:
- `rna_t_ref_count` - RNA reference allele read count
- `rna_t_alt_count` - RNA alternate allele read count
- `rna_t_total_count` - RNA total depth at position
- `rna_t_variant_frequency` - RNA VAF (alt/total)
- `rna_genotyped_variant_frequency` - Same as above (explicitly set)
- `RNA_ID` - RNA sample identifier

Detection logic: check for presence of `rna_t_alt_count` column in the MAF header.

## New Nextflow Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `run_rna_annotation` | `false` | Enable RNA annotation module |
| `run_fusion_neoantigens` | `false` | Enable fusion neoantigen prediction |
| `rna_tpm_threshold` | `1.0` | TPM cutoff for "expressed" classification |

Follows the existing `run_phylowgs` convention for optional pipeline branches.

## Module 1: NEOANTIGEN_RNA_ANNOTATE

### Location

Developed in `mskcc-omics-workflows/modules` repo, installed via nf-core tools.

```
modules/msk/neoantigenutils/rnaannotate/
  main.nf
  meta.yml
  environment.yml
  resources/usr/bin/annotate_rna.py
  tests/
    main.nf.test
    main.nf.test.snap
    nextflow.config
    tags.yml
```

### Inputs

1. Neoantigen TSV (from CONVERTANNOTJSON)
2. Original MAF (same input MAF, may contain rna_* columns)
3. Kallisto abundance.tsv (optional)
4. GTF file (for transcript-to-gene mapping)

### Logic

1. **RNA VAF annotation**: Join neoantigen TSV to MAF on `mutation_id`. The mutation_id format is `chr_pos_ref_alt` (e.g., `7_55259515_T_G`). Reconstruct this key from MAF columns `Chromosome`, `Start_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`. Extract `rna_t_alt_count`, `rna_t_ref_count`, `rna_t_variant_frequency`.

2. **Expression annotation**: Parse Kallisto `abundance.tsv` (columns: `target_id`, `length`, `eff_length`, `est_counts`, `tpm`). Map Ensembl transcript IDs to gene symbols using the GTF file (params.gtf). Annotate each neoantigen with the TPM of its gene. Join key: `Gene` (Hugo_Symbol) in neoantigen TSV to gene symbol from GTF.

3. **Derived boolean**: `rna_expressed` = True if TPM > threshold (configurable) AND (rna_t_alt_count > 0 if VAF data exists).

### Outputs

**Annotated neoantigen TSV** - original columns plus:

| New Column | Source | Description |
|------------|--------|-------------|
| `rna_tpm` | Kallisto | Gene-level TPM value |
| `rna_vaf` | MAF rna_* columns | RNA variant allele frequency |
| `rna_alt_count` | MAF rna_* columns | RNA alt allele read count |
| `rna_ref_count` | MAF rna_* columns | RNA ref allele read count |
| `rna_expressed` | Derived | Boolean: evidence of expression |

**Detailed RNA report TSV** - one row per gene/mutation with full Kallisto and fillout details.

### Data Format Mapping

| FORTE Output | Join Key | Target Column |
|---|---|---|
| MAF `Chromosome:Start_Position:Ref:Alt` | `mutation_id` (chr_pos_ref_alt) | Join key |
| MAF `rna_t_variant_frequency` | - | `rna_vaf` |
| MAF `rna_t_alt_count` | - | `rna_alt_count` |
| MAF `rna_t_ref_count` | - | `rna_ref_count` |
| Kallisto `tpm` (via transcript -> gene) | `Gene` (Hugo_Symbol) | `rna_tpm` |
| Derived | - | `rna_expressed` |

### Mutation ID Reconstruction

The neoantigen pipeline uses a specific mutation_id format constructed in `generate_input.py`:
- SNPs: `{Chromosome}_{Start_Position}_{Reference_Allele}_{Tumor_Seq_Allele2}`
- DELs: `{Chromosome}_{Start_Position}_{Reference_Allele}_D`
- INSs: `{Chromosome}_{Start_Position}_I_{Tumor_Seq_Allele2}`

The annotation script must reconstruct this same format from MAF columns to join correctly.

## Module 2: NEOANTIGEN_FUSION_PREPARE

### Location

```
modules/msk/neoantigenutils/fusionprepare/
  main.nf
  meta.yml
  environment.yml
  resources/usr/bin/prepare_fusion_fasta.py
  tests/
    ...
```

### Inputs

1. AGFusion output directory (per-fusion subdirectories with `.fa` files)
2. HLA string (from existing pipeline)
3. GTF and cDNA files (for NeoSV-compatible formatting)

### Logic

1. Parse AGFusion `.fa` files containing fusion protein sequences
2. Extract peptide windows (9/10/11-mers) around the fusion junction point
3. Create WT equivalents (individual unfused protein segments)
4. Format with FASTA header convention matching the pipeline's encoding scheme
5. Output `*.SV.MUT.fa` and `*.SV.WT.fa` compatible with NETMHCSTABANDPAN input

### Output

Fusion neoantigens appear in a **separate** `fusion_neoantigen.tsv` file, not merged into the main SNV/indel neoantigen TSV, since they have different mutation_id semantics and don't map back to the input MAF.

## Nextflow Channel Architecture

```groovy
// Existing samplesheet (unchanged columns)
ch_samplesheet  // [meta, maf, facets, hla]

// New optional inputs from samplesheet
ch_kallisto     // [meta, kallisto_abundance] - may be empty/null
ch_agfusion     // [meta, agfusion_dir] - may be empty/null

// === EXISTING PIPELINE (unchanged) ===
GENERATE_MUTATED_PEPTIDES(...)
NETMHCSTABANDPAN(...)
NEOANTIGENUTILS_NEOANTIGENINPUT(...)
NEOANTIGEN_EDITING(...)
NEOANTIGENUTILS_CONVERTANNOTJSON(...)

// === NEW: RNA ANNOTATION (optional) ===
if (params.run_rna_annotation) {
    NEOANTIGEN_RNA_ANNOTATE(
        NEOANTIGENUTILS_CONVERTANNOTJSON.out.neoantigenTSV,
        ch_samplesheet.map{ meta, maf, _, _ -> [meta, maf] },
        ch_kallisto
    )
}

// === NEW: FUSION NEOANTIGENS (optional) ===
if (params.run_fusion_neoantigens) {
    NEOANTIGEN_FUSION_PREPARE(ch_agfusion, ch_hla)
    NETMHCSTABANDPAN(fusion_fasta, ...)
    // -> fusion binding predictions -> fusion_neoantigen.tsv
}
```

## Pipeline Output Structure

```
results/
  {sample}/
    neoantigen.tsv                    # Existing (unchanged)
    neoantigen_annotated.json         # Existing (unchanged)
    neoantigen_rna_annotated.tsv      # NEW: main TSV + RNA columns
    rna_detail_report.tsv             # NEW: detailed RNA report
    fusion_neoantigen.tsv             # NEW: fusion neoantigens (separate)
```

## Open Questions for Team Input

1. **Expression threshold**: What TPM cutoff defines "expressed"? Literature suggests TPM > 1 for coding genes. Made configurable via `rna_tpm_threshold` param.
2. **RNA VAF threshold**: What minimum `rna_t_alt_count` constitutes evidence of variant expression? (e.g., >= 2 alt reads?)
3. **Fusion peptide windowing**: AGFusion provides full fusion protein sequences. Window into 9/10/11-mers around the junction, or use full protein context?
4. **Container**: Extend existing `neoantigen-utils-base:1.4.0` or create a new container?
5. **featureCounts vs Kallisto**: FORTE produces both gene-level counts (featureCounts) and transcript-level TPM (Kallisto). Should we support featureCounts as an alternative input?
6. **Kallisto count_features.R**: FORTE also runs `count_features.R` which produces `kallisto.customsummary.txt`. Should we use this processed output or raw `abundance.tsv`?

## References

- Luksza et al. 2017 - Neoantigen quality model
- Luksza et al. 2022 - Updated fitness model
- FORTE pipeline: https://github.com/mskcc/forte (dev branch)
- mskcc-omics-workflows modules: https://github.com/mskcc-omics-workflows/modules
- NeoSV: Structural variant neoantigen prediction
- AGFusion: Annotate gene fusions with protein-level consequences
- GBCMS (GetBaseCountsMultiSample): Variant-level read support from BAMs
