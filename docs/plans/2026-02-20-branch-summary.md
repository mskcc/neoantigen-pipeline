# Branch Summary: `feature/rna-integration-design`

**Date**: 2026-02-20
**Pipeline Repo**: `NeoAntigen_module_update` (branch `feature/rna-integration-design`)
**Module Repo**: `mskcc-omics-workflows/modules` (branch `feature/rna-integration`)

---

## What Was Done

This branch adds optional RNA-seq integration to the neoantigen prediction pipeline. The work spans two repositories and introduces two new capabilities:

1. **RNA Annotation** -- Annotates predicted neoantigens with gene-level expression (TPM from Kallisto) and variant-level RNA support (RNA VAF from FORTE-processed MAFs)
2. **Fusion Neoantigen Preparation** -- Converts AGFusion output into junction-spanning peptide FASTAs suitable for MHC binding prediction

Both features are fully optional. The existing DNA-only pipeline is completely unchanged when these features are disabled (which is the default).

### Commits (Module Repo -- `mskcc-omics-workflows/modules`)

| Commit | Description |
|--------|-------------|
| `601ec6e` | `annotate_rna.py` and unit tests -- core RNA annotation script |
| `ff3fc01` | `prepare_fusion_fasta.py` and unit tests -- fusion peptide extraction |
| `5aae770` | Nextflow modules (`main.nf`, `meta.yml`, `environment.yml`, nf-tests) for both |
| `eb76512` | Fix nf-test snapshots with actual md5 hashes |

### Commits (Pipeline Repo -- `NeoAntigen_module_update`)

| Commit | Description |
|--------|-------------|
| `16704c9` | Design document (`docs/plans/2026-02-18-rna-integration-design.md`) |
| `81396b8` | Implementation plan (`docs/plans/2026-02-18-rna-integration-plan.md`) |
| `af8107a` | Full pipeline integration -- modules, workflow, schema, config (23 files, 1002 insertions) |
| `5f54133` | Test infrastructure -- test data, test config, nf-test file |

### Uncommitted Fixes (Pipeline Repo)

| File | Fix |
|------|-----|
| `annotate_rna.py` | Added `gzip` support for compressed GTF files |
| `tests/rna_annotation.nf.test` | Added `docker` to profile for containerized execution |
| `conf/test_rna.config` | Dynamic samplesheet generation using `${projectDir}` for absolute path resolution |
| `.gitignore` | Ignore dynamically generated `assets/samplesheet_rna_resolved.csv` |
| `tests/rna_annotation.nf.test.snap` | New snapshot file from passing test run |

---

## Architecture: How the RNA Information Is Used

### Design Decision: Annotation, Not Scoring

The neoantigen pipeline predicts neoantigens from DNA mutations and scores them using the Luksza et al. fitness model. That model uses MHC binding affinity (Kd), clonality (logC), and amplitude (logA) to compute a fitness score. **The RNA data does not modify this scoring model.** Instead, RNA information is appended as additional annotation columns for downstream filtering and interpretation.

This was a deliberate architectural choice. The scoring model is validated and published. Changing it would require re-validation. Adding RNA as annotation preserves backward compatibility and lets researchers apply their own RNA-based filters without altering the core predictions.

### Data Flow

```
Existing Pipeline (unchanged):
  MAF + HLA + Facets
    -> MUTALYZER_RETRIEVER
    -> GENERATEMUTFASTA
    -> NETMHCSTABPAN / NETMHC3
    -> FORMATNETMHCPAN
    -> NEOANTIGENINPUT
    -> NEOANTIGEN_EDITING
    -> CONVERTANNOTJSON
    -> neoantigens.tsv (scored predictions with Kd, R, logC, logA, quality)

New Path 1 -- RNA Annotation (runs after CONVERTANNOTJSON):
  neoantigens.tsv + MAF + kallisto_abundance.tsv + GTF
    -> NEOANTIGENUTILS_RNAANNOTATE
    -> neoantigens_rna_annotated.tsv (original columns + rna_tpm, rna_vaf, rna_expressed)
    -> neoantigens_rna_report.tsv   (per-mutation RNA summary)

New Path 2 -- Fusion Neoantigen Preparation:
  agfusion_dir
    -> NEOANTIGENUTILS_FUSIONPREPARE
    -> sample.SV.MUT.fa + sample.SV.WT.fa (junction peptides for MHC binding prediction)
```

### What Each RNA Column Means

The RNA annotation module adds five new columns to the neoantigen TSV:

| Column | Source | Meaning |
|--------|--------|---------|
| `rna_tpm` | Kallisto abundance + GTF | Gene-level expression in transcripts per million. Kallisto quantifies transcript-level TPM; the script maps transcripts to genes via the GTF and sums TPM across all transcripts of a gene. |
| `rna_vaf` | MAF `rna_t_variant_frequency` column | RNA variant allele frequency -- the fraction of RNA reads at the mutation site that carry the variant allele. Provided by FORTE's somatic variant annotation. |
| `rna_alt_count` | MAF `rna_t_alt_count` column | Number of RNA-seq reads supporting the mutant allele at this site. |
| `rna_ref_count` | MAF `rna_t_ref_count` column | Number of RNA-seq reads supporting the reference allele at this site. |
| `rna_expressed` | Computed boolean | `True` if `rna_tpm > threshold` AND (`rna_alt_count` is NA or > 0). This is the key filter column. |

**Graceful degradation**: If Kallisto abundance is not provided, TPM columns are `NA`. If the MAF lacks RNA columns (i.e., the sample was not processed through FORTE's RNA variant calling), VAF columns are `NA`. Both situations result in `rna_expressed = False`.

### How Fusion Peptides Are Extracted

The fusion module parses AGFusion output directories, which contain predicted fusion protein sequences from gene fusion events. For each fusion:

1. Identifies the junction position (marked by `*` in AGFusion output, or estimated)
2. Extracts sliding-window peptides (9, 10, 11 amino acids) spanning the junction
3. Writes paired FASTA files (MUT and WT) compatible with the pipeline's existing NetMHC binding prediction infrastructure

These junction-spanning peptides represent novel amino acid sequences that exist only in the fusion protein and are therefore tumor-specific.

---

## Clinical and Research Utility

### What Problem Does This Solve?

The standard neoantigen prediction pipeline identifies mutations that could produce immunogenic peptides based on DNA sequencing alone. But a DNA mutation that is never transcribed into RNA -- or transcribed at negligible levels -- cannot produce a protein that the immune system can recognize. Without RNA data, the pipeline has no way to distinguish between:

- A mutation in a highly expressed gene (likely to produce targetable peptides)
- A mutation in a silenced gene (will never produce a targetable peptide)
- A mutation where the variant allele is present in DNA but absent from RNA (allelic imbalance, monoallelic expression, or nonsense-mediated decay)

### Value to Researchers

**Neoantigen prioritization.** The pipeline typically predicts hundreds of candidate neoantigens. Researchers need to prioritize a small number for experimental validation (e.g., T-cell reactivity assays, tetramer staining, ELISpot). RNA expression is one of the strongest filters available -- published work shows that integrating RNA-seq data reduces candidates while enriching for immunogenic neoantigens ([Frontiers in Immunology, 2023](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1251603/full); [Cancer Research, PMC5329159](https://pmc.ncbi.nlm.nih.gov/articles/PMC5329159/)).

**Fusion neoantigens are an underexplored source.** Gene fusions drive ~16% of cancers and produce novel junction peptides that are absent from normal tissue. Fusion-derived neoantigens tend to be more immunogenic, have more targets per event, and are more likely to be shared across patients compared to SNV-derived neoantigens ([Yang et al., Cancer Letters, 2021](https://pubmed.ncbi.nlm.nih.gov/33675984/); [Mertens et al., Nature Medicine, 2019](https://www.nature.com/articles/s41591-019-0434-2)). By extracting junction peptides and routing them through existing MHC binding prediction, this pipeline enables researchers to evaluate fusion neoantigens alongside point-mutation neoantigens in a unified framework.

**Multi-omic evidence integration.** Having DNA-predicted neoantigens annotated with RNA evidence in a single TSV enables researchers to apply composite filters (e.g., "high MHC binding affinity AND expressed AND variant detected in RNA") without writing custom join scripts or switching between tools.

### Value to Clinicians

**Personalized vaccine design.** For neoantigen-based cancer vaccines (now in clinical trials at multiple centers), clinicians must select a small number of peptides (typically 10-20) from often hundreds of candidates. RNA expression data directly answers: "Is this mutation actually producing protein in the patient's tumor?" A neoantigen vaccine targeting a silenced gene is wasted effort.

**Immunotherapy response prediction.** Tumor mutational burden (TMB) is an imperfect predictor of checkpoint inhibitor response. The detection of mutant mRNA in the tumor is increasingly recognized as a better biomarker than raw neoantigen counts. A VAF cutoff of >= 0.04 in RNA-seq data has been shown to yield high positive predictive rates for true mutant expression ([Frontiers in Immunology, 2023](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1301100/full)).

**Fusion-targeted therapy.** Some fusion neoantigens are shared across patients (e.g., DNAJB1-PRKACA in fibrolamellar carcinoma, SYT-SSX in synovial sarcoma). Phase I trials targeting fusion neoantigens have shown tumor-specific T-cell responses ([Cell Reports Medicine, 2024](https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(24)00115-0); [Frontiers in Oncology, 2024](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2024.1367450/full)).

### Value to Patients

The ultimate impact is better-targeted therapy. By filtering out neoantigens that are not expressed, clinicians can focus on candidates more likely to elicit an immune response. This means:

- **More effective personalized vaccines** -- peptides selected based on multi-omic evidence rather than DNA alone
- **Fewer wasted therapy slots** -- neoantigen vaccines typically carry a limited number of peptides; excluding silent mutations makes room for more promising targets
- **Access to fusion-derived targets** -- patients with fusion-driven cancers (sarcomas, certain carcinomas, pediatric cancers) gain access to a neoantigen source that was previously not evaluated by this pipeline

---

## Is the Data Actually Useful? A Candid Assessment

### What the evidence supports

**RNA expression filtering is well-supported.** Multiple studies demonstrate that neoantigens in unexpressed genes do not elicit T-cell responses, and that RNA-based filtering improves the hit rate of validation experiments. This is one of the least controversial filters in the neoantigen prediction field. The TPM threshold approach used here (default 1.0) is standard practice.

**RNA VAF adds independent evidence.** Even if a gene is expressed, the mutant allele may not be transcribed (e.g., the mutation is on the silenced allele, or nonsense-mediated decay eliminates the mutant transcript). RNA VAF directly addresses this. However, RNA VAF is noisy -- low-coverage sites can produce unreliable estimates, and the absence of RNA reads at a site does not definitively prove the allele is not expressed.

**Fusion neoantigens are genuinely promising but early.** The immunogenicity advantage of fusion neoantigens is supported by Nature Medicine-published data. However, the field is still young -- most evidence comes from case reports and small trials. The module correctly focuses on junction peptides (the truly novel sequences) rather than full fusion proteins.

### Limitations of this implementation

1. **Annotation only, not scoring.** The `rna_expressed` boolean is a simple threshold. It does not integrate into the Luksza fitness model. Researchers still need to apply their own prioritization logic combining RNA evidence with binding affinity and clonality. This is intentional (preserves the validated model) but means the RNA data requires interpretation.

2. **Gene-level TPM, not allele-specific.** Kallisto provides total gene expression. If a gene has high TPM but the mutation is on the unexpressed allele, TPM alone would be misleading. The RNA VAF column partially addresses this, but only when the MAF contains FORTE RNA columns.

3. **No RNA-based binding affinity change.** The pipeline does not re-score binding affinity based on RNA expression levels. Expression and binding are treated as independent axes.

4. **Fusion MUT and WT are identical.** The current fusion module writes the same junction peptide as both MUT and WT. This is because the "wildtype" for a fusion junction does not have a natural analog -- the junction sequence simply does not exist in normal cells. This may affect downstream tools that compare MUT vs WT binding.

5. **Depends on upstream FORTE processing.** RNA annotation quality is entirely dependent on the quality of the FORTE pipeline outputs. Bad Kallisto quantification or incorrect RNA variant calling will propagate.

### Bottom line

The RNA annotation adds genuinely useful information that the DNA-only pipeline cannot provide. Expression filtering is the single most impactful filter for reducing false-positive neoantigens, and this implementation follows established practice. The fusion neoantigen preparation opens up an entirely new class of targets. Neither feature changes the existing pipeline's behavior unless explicitly enabled, so there is no risk to current users.

The data is useful. It is not sufficient on its own for clinical decisions -- it is one layer in a multi-layered prioritization process -- but it fills a real gap in the current pipeline's output.

---

## Files Changed (Complete List)

### Module Repo (`mskcc-omics-workflows/modules`, branch `feature/rna-integration`)

**New files:**
- `modules/msk/neoantigenutils/rnaannotate/main.nf` -- Nextflow process definition
- `modules/msk/neoantigenutils/rnaannotate/meta.yml` -- Module metadata
- `modules/msk/neoantigenutils/rnaannotate/environment.yml` -- Conda environment
- `modules/msk/neoantigenutils/rnaannotate/resources/usr/bin/annotate_rna.py` -- Core annotation script
- `modules/msk/neoantigenutils/rnaannotate/resources/usr/bin/test_annotate_rna.py` -- Unit tests
- `modules/msk/neoantigenutils/rnaannotate/tests/main.nf.test` -- nf-test
- `modules/msk/neoantigenutils/rnaannotate/tests/main.nf.test.snap` -- nf-test snapshot
- `modules/msk/neoantigenutils/rnaannotate/tests/tags.yml` -- nf-test tags
- `modules/msk/neoantigenutils/fusionprepare/main.nf` -- Nextflow process definition
- `modules/msk/neoantigenutils/fusionprepare/meta.yml` -- Module metadata
- `modules/msk/neoantigenutils/fusionprepare/environment.yml` -- Conda environment
- `modules/msk/neoantigenutils/fusionprepare/resources/usr/bin/prepare_fusion_fasta.py` -- Fusion script
- `modules/msk/neoantigenutils/fusionprepare/resources/usr/bin/test_prepare_fusion.py` -- Unit tests
- `modules/msk/neoantigenutils/fusionprepare/tests/main.nf.test` -- nf-test
- `modules/msk/neoantigenutils/fusionprepare/tests/main.nf.test.snap` -- nf-test snapshot
- `modules/msk/neoantigenutils/fusionprepare/tests/tags.yml` -- nf-test tags

### Pipeline Repo (`NeoAntigen_module_update`, branch `feature/rna-integration-design`)

**New files:**
- `docs/plans/2026-02-18-rna-integration-design.md` -- Design document
- `docs/plans/2026-02-18-rna-integration-plan.md` -- Implementation plan
- `modules/msk/neoantigenutils/rnaannotate/` -- Entire module (copied from module repo)
- `modules/msk/neoantigenutils/fusionprepare/` -- Entire module (copied from module repo)
- `tests/data/test_abundance.tsv` -- Fake Kallisto abundance for testing
- `tests/rna_annotation.nf.test` -- Pipeline-level nf-test
- `tests/rna_annotation.nf.test.snap` -- Passing snapshot
- `assets/samplesheet_rna.csv` -- RNA test samplesheet template
- `conf/test_rna.config` -- Test profile configuration

**Modified files:**
- `workflows/neoantigenpipeline.nf` -- Added RNA/fusion channels, conditional processing blocks, imports
- `assets/schema_input.json` -- Added optional `kallisto_abundance` and `agfusion_dir` columns
- `nextflow_schema.json` -- Added `run_rna_annotation`, `run_fusion_neoantigens`, `rna_tpm_threshold` params
- `nextflow.config` -- Default param values + `test_rna` profile
- `conf/modules.config` -- RNAANNOTATE process config with TPM threshold passthrough
- `modules.json` -- Module registry entries for both new modules
- `subworkflows/local/utils_mskcc_neoantigenpipeline_pipeline/main.nf` -- Samplesheet parsing for 6 columns
- `nf-test.config` -- Added `conf/test_rna.config` to triggers
- `tests/.nftignore` -- Added RNA report pattern
- `tests/nextflow.config` -- Added RNA test data reference
- `.gitignore` -- Added exceptions for test data and generated samplesheet
