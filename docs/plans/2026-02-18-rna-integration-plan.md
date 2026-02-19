# RNA-seq Integration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add optional RNA annotation (expression TPM, RNA VAF) and fusion neoantigen prediction to the neoantigen pipeline, consuming pre-computed FORTE outputs.

**Architecture:** Post-processing annotation module runs after existing pipeline. RNA data adds columns to output TSV without changing the scoring model. Fusion neoantigens go through existing MHC binding machinery into a separate output file.

**Tech Stack:** Nextflow DSL2, Python 3, pandas, pyensembl, nf-core module conventions

**Design doc:** `docs/plans/2026-02-18-rna-integration-design.md`

---

## Two Repositories

This implementation spans two repos:

1. **Module repo** (`mskcc-omics-workflows/modules`) — where new Nextflow modules are developed and tested
2. **Pipeline repo** (`NeoAntigen_module_update`) — where modules are installed and integrated into the workflow

The plan covers both. Module development (Tasks 1-6) happens in the module repo. Pipeline integration (Tasks 7-12) happens in this repo after modules are published.

---

### Task 1: Create `annotate_rna.py` — Mutation ID Builder + MAF RNA Column Extraction

**Context:** This is the core annotation script. It must reconstruct the neoantigen pipeline's mutation_id format from MAF columns, then join RNA columns from the MAF onto the neoantigen TSV.

**Files:**
- Create: `modules/msk/neoantigenutils/rnaannotate/resources/usr/bin/annotate_rna.py`
- Create: `modules/msk/neoantigenutils/rnaannotate/resources/usr/bin/test_annotate_rna.py`

**Step 1: Write test for mutation_id reconstruction**

The mutation_id format (from `generate_input.py` lines 132-274):
- SNP: `{Chromosome}_{Start_Position}_{Reference_Allele}_{Tumor_Seq_Allele2}`
- DEL: `{Chromosome}_{Start_Position}_{Reference_Allele}_D`
- INS: `{Chromosome}_{Start_Position}_I_{Tumor_Seq_Allele2}`

```python
# test_annotate_rna.py
import pytest
import pandas as pd
import tempfile
import os
from annotate_rna import build_mutation_id, extract_rna_columns_from_maf, annotate_expression, merge_annotations

class TestBuildMutationId:
    def test_snp(self):
        row = {"Chromosome": "7", "Start_Position": 55259515,
               "Reference_Allele": "T", "Tumor_Seq_Allele2": "G",
               "Variant_Type": "SNP"}
        assert build_mutation_id(row) == "7_55259515_T_G"

    def test_del(self):
        row = {"Chromosome": "17", "Start_Position": 7577120,
               "Reference_Allele": "AC", "Tumor_Seq_Allele2": "-",
               "Variant_Type": "DEL"}
        assert build_mutation_id(row) == "17_7577120_AC_D"

    def test_ins(self):
        row = {"Chromosome": "12", "Start_Position": 25398284,
               "Reference_Allele": "-", "Tumor_Seq_Allele2": "AGT",
               "Variant_Type": "INS"}
        assert build_mutation_id(row) == "12_25398284_I_AGT"

    def test_dnp(self):
        row = {"Chromosome": "1", "Start_Position": 115256530,
               "Reference_Allele": "TT", "Tumor_Seq_Allele2": "AA",
               "Variant_Type": "DNP"}
        assert build_mutation_id(row) == "1_115256530_TT_AA"


class TestExtractRnaColumns:
    def test_maf_with_rna_columns(self):
        maf_df = pd.DataFrame({
            "Chromosome": ["7", "17"],
            "Start_Position": [55259515, 7577120],
            "Reference_Allele": ["T", "AC"],
            "Tumor_Seq_Allele2": ["G", "-"],
            "Variant_Type": ["SNP", "DEL"],
            "Hugo_Symbol": ["EGFR", "TP53"],
            "rna_t_alt_count": [15, 0],
            "rna_t_ref_count": [85, 120],
            "rna_t_variant_frequency": [0.15, 0.0],
        })
        result = extract_rna_columns_from_maf(maf_df)
        assert len(result) == 2
        assert result.loc[result["mutation_id"] == "7_55259515_T_G", "rna_alt_count"].values[0] == 15
        assert result.loc[result["mutation_id"] == "17_7577120_AC_D", "rna_vaf"].values[0] == 0.0

    def test_maf_without_rna_columns(self):
        maf_df = pd.DataFrame({
            "Chromosome": ["7"],
            "Start_Position": [55259515],
            "Reference_Allele": ["T"],
            "Tumor_Seq_Allele2": ["G"],
            "Variant_Type": ["SNP"],
            "Hugo_Symbol": ["EGFR"],
        })
        result = extract_rna_columns_from_maf(maf_df)
        assert len(result) == 1
        assert pd.isna(result["rna_alt_count"].values[0])
        assert pd.isna(result["rna_vaf"].values[0])


class TestAnnotateExpression:
    def test_kallisto_to_gene_tpm(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("target_id\tlength\teff_length\test_counts\ttpm\n")
            f.write("ENST00000275493.7\t3188\t2939.42\t1500\t45.2\n")
            f.write("ENST00000454757.1\t637\t388.42\t50\t1.1\n")
            abundance_path = f.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.gtf', delete=False) as f:
            f.write('chr7\tensembl\ttranscript\t55019017\t55211628\t.\t+\t.\tgene_id "ENSG00000146648"; transcript_id "ENST00000275493"; gene_name "EGFR";\n')
            f.write('chr7\tensembl\ttranscript\t55019017\t55211628\t.\t+\t.\tgene_id "ENSG00000146648"; transcript_id "ENST00000454757"; gene_name "EGFR";\n')
            gtf_path = f.name

        result = annotate_expression(abundance_path, gtf_path)
        # EGFR should have summed TPM from both transcripts
        assert "EGFR" in result["Gene"].values
        egfr_tpm = result.loc[result["Gene"] == "EGFR", "rna_tpm"].values[0]
        assert abs(egfr_tpm - 46.3) < 0.1  # 45.2 + 1.1

        os.unlink(abundance_path)
        os.unlink(gtf_path)

    def test_no_kallisto_returns_empty(self):
        result = annotate_expression(None, None)
        assert result is None


class TestMergeAnnotations:
    def test_full_merge(self):
        neo_tsv = pd.DataFrame({
            "id": ["neo1", "neo2"],
            "mutation_id": ["7_55259515_T_G", "17_7577120_AC_D"],
            "Gene": ["EGFR", "TP53"],
            "quality": [0.85, 0.42],
        })
        rna_vaf = pd.DataFrame({
            "mutation_id": ["7_55259515_T_G", "17_7577120_AC_D"],
            "rna_alt_count": [15, 0],
            "rna_ref_count": [85, 120],
            "rna_vaf": [0.15, 0.0],
        })
        expression = pd.DataFrame({
            "Gene": ["EGFR", "TP53"],
            "rna_tpm": [46.3, 0.2],
        })
        result = merge_annotations(neo_tsv, rna_vaf, expression, tpm_threshold=1.0)
        assert "rna_expressed" in result.columns
        # EGFR: TPM > 1 AND rna_alt_count > 0 => True
        assert result.loc[result["Gene"] == "EGFR", "rna_expressed"].values[0] == True
        # TP53: TPM < 1 => False
        assert result.loc[result["Gene"] == "TP53", "rna_expressed"].values[0] == False
```

**Step 2: Run tests to verify they fail**

Run: `cd modules/msk/neoantigenutils/rnaannotate/resources/usr/bin && python -m pytest test_annotate_rna.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'annotate_rna'`

**Step 3: Implement `annotate_rna.py`**

```python
#!/usr/bin/env python3
"""Annotate neoantigen TSV with RNA-seq data from FORTE outputs."""

import argparse
import re
import sys

import pandas as pd

VERSION = "1.0.0"


def build_mutation_id(row):
    """Reconstruct mutation_id from MAF columns.

    Matches the format used by generate_input.py:
    - SNP/DNP/TNP: chr_pos_ref_alt
    - DEL: chr_pos_ref_D
    - INS: chr_pos_I_alt
    """
    chrom = str(row["Chromosome"])
    pos = str(row["Start_Position"])
    ref = str(row["Reference_Allele"])
    alt = str(row["Tumor_Seq_Allele2"])
    vtype = str(row["Variant_Type"])

    if vtype == "DEL":
        return f"{chrom}_{pos}_{ref}_D"
    elif vtype == "INS":
        return f"{chrom}_{pos}_I_{alt}"
    else:
        return f"{chrom}_{pos}_{ref}_{alt}"


def extract_rna_columns_from_maf(maf_df):
    """Extract RNA VAF columns from MAF if present.

    Returns DataFrame with mutation_id, rna_alt_count, rna_ref_count, rna_vaf.
    If rna columns are absent, returns DataFrame with NaN values.
    """
    maf_df = maf_df.copy()
    maf_df["mutation_id"] = maf_df.apply(build_mutation_id, axis=1)

    has_rna = "rna_t_alt_count" in maf_df.columns

    result = pd.DataFrame({"mutation_id": maf_df["mutation_id"]})

    if has_rna:
        result["rna_alt_count"] = maf_df["rna_t_alt_count"].values
        result["rna_ref_count"] = maf_df["rna_t_ref_count"].values
        result["rna_vaf"] = maf_df["rna_t_variant_frequency"].values
    else:
        result["rna_alt_count"] = pd.NA
        result["rna_ref_count"] = pd.NA
        result["rna_vaf"] = pd.NA

    return result


def parse_gtf_gene_map(gtf_path):
    """Parse GTF to build transcript_id -> gene_name mapping."""
    tx2gene = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8]
            tx_match = re.search(r'transcript_id "([^"]+)"', attrs)
            gene_match = re.search(r'gene_name "([^"]+)"', attrs)
            if tx_match and gene_match:
                tx_id = tx_match.group(1).split(".")[0]
                gene_name = gene_match.group(1)
                tx2gene[tx_id] = gene_name
    return tx2gene


def annotate_expression(abundance_path, gtf_path):
    """Parse Kallisto abundance and map to gene-level TPM.

    Returns DataFrame with Gene, rna_tpm columns.
    Returns None if abundance_path is None.
    """
    if abundance_path is None:
        return None

    abundance = pd.read_csv(abundance_path, sep="\t")
    tx2gene = parse_gtf_gene_map(gtf_path)

    # Strip version from target_id
    abundance["transcript_id"] = abundance["target_id"].str.split(".").str[0]
    abundance["Gene"] = abundance["transcript_id"].map(tx2gene)

    # Sum TPM per gene
    gene_tpm = (
        abundance.dropna(subset=["Gene"])
        .groupby("Gene")["tpm"]
        .sum()
        .reset_index()
        .rename(columns={"tpm": "rna_tpm"})
    )

    return gene_tpm


def merge_annotations(neo_tsv, rna_vaf_df, expression_df, tpm_threshold=1.0):
    """Merge RNA annotations onto neoantigen TSV.

    Returns annotated DataFrame with rna_tpm, rna_vaf, rna_alt_count,
    rna_ref_count, rna_expressed columns.
    """
    result = neo_tsv.copy()

    # Merge RNA VAF
    if rna_vaf_df is not None:
        result = result.merge(rna_vaf_df, on="mutation_id", how="left")
    else:
        result["rna_alt_count"] = pd.NA
        result["rna_ref_count"] = pd.NA
        result["rna_vaf"] = pd.NA

    # Merge expression
    if expression_df is not None:
        result = result.merge(expression_df, on="Gene", how="left")
    else:
        result["rna_tpm"] = pd.NA

    # Compute rna_expressed
    tpm_ok = result["rna_tpm"].notna() & (result["rna_tpm"] > tpm_threshold)
    vaf_ok = result["rna_alt_count"].isna() | (result["rna_alt_count"] > 0)
    result["rna_expressed"] = tpm_ok & vaf_ok

    return result


def main():
    parser = argparse.ArgumentParser(description="Annotate neoantigens with RNA data")
    parser.add_argument("--neoantigen_tsv", required=True, help="Neoantigen TSV from convertannotjson")
    parser.add_argument("--maf", required=True, help="Input MAF (may contain rna_* columns)")
    parser.add_argument("--kallisto_abundance", default=None, help="Kallisto abundance.tsv (optional)")
    parser.add_argument("--gtf", default=None, help="GTF file for transcript-to-gene mapping")
    parser.add_argument("--tpm_threshold", type=float, default=1.0, help="TPM threshold for expressed")
    parser.add_argument("--output_annotated", required=True, help="Output annotated TSV")
    parser.add_argument("--output_report", required=True, help="Output detailed RNA report")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()

    neo_tsv = pd.read_csv(args.neoantigen_tsv, sep="\t")
    maf_df = pd.read_csv(args.maf, sep="\t", comment="#")

    rna_vaf_df = extract_rna_columns_from_maf(maf_df)
    expression_df = annotate_expression(args.kallisto_abundance, args.gtf)

    annotated = merge_annotations(neo_tsv, rna_vaf_df, expression_df, args.tpm_threshold)
    annotated.to_csv(args.output_annotated, sep="\t", index=False)

    # Detailed report: one row per gene with all RNA details
    report_cols = ["mutation_id", "Gene", "rna_tpm", "rna_vaf", "rna_alt_count", "rna_ref_count", "rna_expressed"]
    available_cols = [c for c in report_cols if c in annotated.columns]
    report = annotated[available_cols].drop_duplicates()
    report.to_csv(args.output_report, sep="\t", index=False)


if __name__ == "__main__":
    main()
```

**Step 4: Run tests to verify they pass**

Run: `cd modules/msk/neoantigenutils/rnaannotate/resources/usr/bin && python -m pytest test_annotate_rna.py -v`
Expected: All PASS

**Step 5: Commit**

```bash
git add modules/msk/neoantigenutils/rnaannotate/resources/usr/bin/annotate_rna.py \
       modules/msk/neoantigenutils/rnaannotate/resources/usr/bin/test_annotate_rna.py
git commit -m "feat: add annotate_rna.py with tests for RNA annotation"
```

---

### Task 2: Create Nextflow Module for RNA Annotation

**Context:** Standard nf-core module structure following existing conventions in the neoantigen pipeline.

**Files:**
- Create: `modules/msk/neoantigenutils/rnaannotate/main.nf`
- Create: `modules/msk/neoantigenutils/rnaannotate/meta.yml`
- Create: `modules/msk/neoantigenutils/rnaannotate/environment.yml`

**Step 1: Write `main.nf`**

Model after `modules/msk/neoantigenutils/convertannotjson/main.nf` (same container, similar structure):

```groovy
process NEOANTIGENUTILS_RNAANNOTATE {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/mskcc-omics-workflows/neoantigen-utils-base:1.4.0':
        'ghcr.io/mskcc-omics-workflows/neoantigen-utils-base:1.4.0' }"

    input:
    tuple val(meta), path(neoantigen_tsv)
    tuple val(meta2), path(maf)
    tuple val(meta3), path(kallisto_abundance)
    path(gtf)

    output:
    tuple val(meta), path("*_rna_annotated.tsv"),  emit: annotated_tsv
    tuple val(meta), path("*_rna_report.tsv"),     emit: rna_report
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kallisto_arg = kallisto_abundance ? "--kallisto_abundance ${kallisto_abundance} --gtf ${gtf}" : ""

    """
    annotate_rna.py \\
        --neoantigen_tsv ${neoantigen_tsv} \\
        --maf ${maf} \\
        ${kallisto_arg} \\
        --output_annotated ${prefix}_rna_annotated.tsv \\
        --output_report ${prefix}_rna_report.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotate_rna: \$(annotate_rna.py -v | sed 's/annotate_rna.py //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rna_annotated.tsv
    touch ${prefix}_rna_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotate_rna: ${VERSION}
    END_VERSIONS
    """
}
```

**Step 2: Write `meta.yml`**

```yaml
name: "neoantigenutils_rnaannotate"
description: Annotate neoantigen TSV with RNA expression and VAF data from FORTE outputs
keywords:
  - neoantigen
  - rna
  - expression
  - annotation
tools:
  - annotate_rna:
      description: Annotates neoantigens with RNA TPM and VAF from FORTE pipeline outputs
      homepage: https://github.com/mskcc-omics-workflows/modules
      licence: ["MIT"]

input:
  - meta:
      type: map
      description: Sample metadata
  - neoantigen_tsv:
      type: file
      description: Neoantigen TSV from convertannotjson
      pattern: "*.tsv"
  - maf:
      type: file
      description: Input MAF file (may contain rna_* columns from FORTE fillout)
      pattern: "*.maf"
  - kallisto_abundance:
      type: file
      description: Kallisto abundance.tsv with transcript-level TPM (optional)
      pattern: "abundance.tsv"
  - gtf:
      type: file
      description: GTF annotation file for transcript-to-gene mapping
      pattern: "*.gtf"

output:
  - annotated_tsv:
      type: file
      description: Neoantigen TSV with RNA annotation columns added
      pattern: "*_rna_annotated.tsv"
  - rna_report:
      type: file
      description: Detailed RNA report per gene/mutation
      pattern: "*_rna_report.tsv"
  - versions:
      type: file
      description: File containing software versions
      pattern: "versions.yml"
```

**Step 3: Write `environment.yml`**

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.10
  - pandas=2.0
```

**Step 4: Commit**

```bash
git add modules/msk/neoantigenutils/rnaannotate/main.nf \
       modules/msk/neoantigenutils/rnaannotate/meta.yml \
       modules/msk/neoantigenutils/rnaannotate/environment.yml
git commit -m "feat: add Nextflow module definition for RNA annotation"
```

---

### Task 3: Create `prepare_fusion_fasta.py` with Tests

**Context:** This script converts AGFusion output (fusion protein FASTAs) into the format expected by the neoantigen pipeline's NETMHCSTABANDPAN module. It extracts peptide windows around the fusion junction.

**Files:**
- Create: `modules/msk/neoantigenutils/fusionprepare/resources/usr/bin/prepare_fusion_fasta.py`
- Create: `modules/msk/neoantigenutils/fusionprepare/resources/usr/bin/test_prepare_fusion.py`

**Step 1: Write tests**

```python
# test_prepare_fusion.py
import pytest
import tempfile
import os
from prepare_fusion_fasta import parse_agfusion_dir, extract_junction_peptides, write_fasta_pair

class TestExtractJunctionPeptides:
    def test_9mer_window(self):
        # Fusion protein: ...ABCDEFGH|IJKLMNOP...  (| = junction at position 8)
        fusion_seq = "ABCDEFGHIJKLMNOP"
        junction_pos = 8
        peptides = extract_junction_peptides(fusion_seq, junction_pos, peptide_lengths=[9])
        # 9-mers spanning junction: must include at least 1 AA from each side
        assert len(peptides) > 0
        for pep in peptides:
            assert len(pep) == 9

    def test_multiple_lengths(self):
        fusion_seq = "ABCDEFGHIJKLMNOPQRSTUVWX"
        junction_pos = 12
        peptides = extract_junction_peptides(fusion_seq, junction_pos, peptide_lengths=[9, 10, 11])
        lengths = set(len(p) for p in peptides)
        assert 9 in lengths
        assert 10 in lengths
        assert 11 in lengths

    def test_short_sequence_handled(self):
        fusion_seq = "ABCD"
        junction_pos = 2
        peptides = extract_junction_peptides(fusion_seq, junction_pos, peptide_lengths=[9])
        # Too short for 9-mer, should return empty
        assert len(peptides) == 0


class TestParseAgfusionDir:
    def test_reads_fa_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fusion_dir = os.path.join(tmpdir, "EGFR-SEPT14")
            os.makedirs(fusion_dir)
            with open(os.path.join(fusion_dir, "EGFR-SEPT14_protein.fa"), "w") as f:
                f.write(">ENST00000275493-ENST00000400454\n")
                f.write("MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRD\n")

            fusions = parse_agfusion_dir(tmpdir)
            assert len(fusions) >= 1
            assert fusions[0]["gene_pair"] == "EGFR-SEPT14"


class TestWriteFastaPair:
    def test_creates_mut_and_wt_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            peptides = [
                {"id": "EGFR_SEPT14_1", "mut_seq": "ABCDEFGHI", "wt_seq": "ABCXYZGHI"}
            ]
            mut_path = os.path.join(tmpdir, "test.SV.MUT.fa")
            wt_path = os.path.join(tmpdir, "test.SV.WT.fa")
            write_fasta_pair(peptides, mut_path, wt_path)
            assert os.path.exists(mut_path)
            assert os.path.exists(wt_path)
            with open(mut_path) as f:
                content = f.read()
                assert ">EGFR_SEPT14_1_M" in content
                assert "ABCDEFGHI" in content
```

**Step 2: Run tests — expect FAIL**

**Step 3: Implement `prepare_fusion_fasta.py`**

```python
#!/usr/bin/env python3
"""Convert AGFusion output to pipeline-compatible FASTA format."""

import argparse
import os
import re

VERSION = "1.0.0"


def parse_agfusion_dir(agfusion_dir):
    """Parse AGFusion output directory for fusion protein sequences.

    AGFusion creates subdirectories per fusion event, each containing
    *_protein.fa files with fusion protein sequences.
    """
    fusions = []
    for entry in os.listdir(agfusion_dir):
        subdir = os.path.join(agfusion_dir, entry)
        if not os.path.isdir(subdir):
            continue
        gene_pair = entry  # e.g., "EGFR-SEPT14"
        for fname in os.listdir(subdir):
            if fname.endswith("_protein.fa"):
                fpath = os.path.join(subdir, fname)
                sequences = read_fasta(fpath)
                for seq_id, seq in sequences:
                    fusions.append({
                        "gene_pair": gene_pair,
                        "transcript_pair": seq_id,
                        "sequence": seq,
                    })
    return fusions


def read_fasta(path):
    """Read FASTA file, return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line)
    if current_header:
        sequences.append((current_header, "".join(current_seq)))
    return sequences


def extract_junction_peptides(fusion_seq, junction_pos, peptide_lengths=None):
    """Extract peptide windows spanning the fusion junction.

    For each peptide length, slide a window across the junction point.
    Each peptide must include at least 1 AA from each fusion partner.
    """
    if peptide_lengths is None:
        peptide_lengths = [9, 10, 11]
    peptides = []
    for plen in peptide_lengths:
        if len(fusion_seq) < plen:
            continue
        # Window must span junction: start from max(0, junction-plen+1) to min(junction, len-plen)
        start_min = max(0, junction_pos - plen + 1)
        start_max = min(junction_pos, len(fusion_seq) - plen)
        for start in range(start_min, start_max + 1):
            pep = fusion_seq[start:start + plen]
            if len(pep) == plen:
                peptides.append(pep)
    return peptides


def infer_junction_position(fusion_seq, five_prime_len=None):
    """Infer junction position from AGFusion sequence.

    AGFusion marks the junction with '*' in some outputs.
    If not present, use five_prime_len if provided, else midpoint.
    """
    if "*" in fusion_seq:
        return fusion_seq.index("*")
    if five_prime_len is not None:
        return five_prime_len
    return len(fusion_seq) // 2


def write_fasta_pair(peptides, mut_path, wt_path):
    """Write MUT and WT FASTA files for pipeline compatibility."""
    with open(mut_path, "w") as mut_f, open(wt_path, "w") as wt_f:
        for pep in peptides:
            mut_f.write(f">{pep['id']}_M\n{pep['mut_seq']}\n")
            wt_f.write(f">{pep['id']}_W\n{pep['wt_seq']}\n")


def main():
    parser = argparse.ArgumentParser(description="Prepare fusion FASTAs for neoantigen pipeline")
    parser.add_argument("--agfusion_dir", required=True, help="AGFusion output directory")
    parser.add_argument("--output_prefix", required=True, help="Output file prefix")
    parser.add_argument("--peptide_lengths", default="9,10,11", help="Comma-separated peptide lengths")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()

    peptide_lengths = [int(x) for x in args.peptide_lengths.split(",")]
    fusions = parse_agfusion_dir(args.agfusion_dir)

    all_peptides = []
    pep_counter = 0
    for fusion in fusions:
        seq = fusion["sequence"].replace("*", "")
        junction = infer_junction_position(fusion["sequence"])
        junction_peps = extract_junction_peptides(seq, junction, peptide_lengths)
        for pep in junction_peps:
            pep_counter += 1
            all_peptides.append({
                "id": f"{fusion['gene_pair']}_{pep_counter}",
                "mut_seq": pep,
                "wt_seq": pep,  # Placeholder — WT determination requires additional logic
            })

    mut_path = f"{args.output_prefix}.SV.MUT.fa"
    wt_path = f"{args.output_prefix}.SV.WT.fa"
    write_fasta_pair(all_peptides, mut_path, wt_path)


if __name__ == "__main__":
    main()
```

**Step 4: Run tests — expect PASS**

**Step 5: Commit**

```bash
git add modules/msk/neoantigenutils/fusionprepare/resources/usr/bin/prepare_fusion_fasta.py \
       modules/msk/neoantigenutils/fusionprepare/resources/usr/bin/test_prepare_fusion.py
git commit -m "feat: add prepare_fusion_fasta.py with tests"
```

---

### Task 4: Create Nextflow Module for Fusion Preparation

**Files:**
- Create: `modules/msk/neoantigenutils/fusionprepare/main.nf`
- Create: `modules/msk/neoantigenutils/fusionprepare/meta.yml`
- Create: `modules/msk/neoantigenutils/fusionprepare/environment.yml`

**Step 1: Write `main.nf`**

```groovy
process NEOANTIGENUTILS_FUSIONPREPARE {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/mskcc-omics-workflows/neoantigen-utils-base:1.4.0':
        'ghcr.io/mskcc-omics-workflows/neoantigen-utils-base:1.4.0' }"

    input:
    tuple val(meta), path(agfusion_dir)

    output:
    tuple val(meta), path("*.SV.MUT.fa"), emit: mut_fasta
    tuple val(meta), path("*.SV.WT.fa"),  emit: wt_fasta
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prepare_fusion_fasta.py \\
        --agfusion_dir ${agfusion_dir} \\
        --output_prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare_fusion_fasta: \$(prepare_fusion_fasta.py -v | sed 's/prepare_fusion_fasta.py //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.SV.MUT.fa
    touch ${prefix}.SV.WT.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prepare_fusion_fasta: ${VERSION}
    END_VERSIONS
    """
}
```

**Step 2: Write `meta.yml` and `environment.yml`** (similar to Task 2, adapted for fusion)

**Step 3: Commit**

```bash
git add modules/msk/neoantigenutils/fusionprepare/
git commit -m "feat: add Nextflow module for fusion FASTA preparation"
```

---

### Task 5: Update Samplesheet Schema

**Context:** Add optional `kallisto_abundance` and `agfusion_dir` columns to the samplesheet schema.

**Files:**
- Modify: `assets/schema_input.json`

**Step 1: Add new optional properties to schema**

Add after the `hla_file` property (line 36):

```json
"kallisto_abundance": {
    "type": "string",
    "format": "file-path",
    "pattern": "^\\S+\\.tsv$",
    "errorMessage": "Kallisto abundance.tsv from FORTE pipeline (optional)"
},
"agfusion_dir": {
    "type": "string",
    "format": "directory-path",
    "errorMessage": "AGFusion output directory from FORTE pipeline (optional)"
}
```

Note: These are NOT added to the `required` array (line 38).

**Step 2: Commit**

```bash
git add assets/schema_input.json
git commit -m "feat: add optional RNA columns to samplesheet schema"
```

---

### Task 6: Update Pipeline Parameters Schema

**Context:** Add `run_rna_annotation`, `run_fusion_neoantigens`, and `rna_tpm_threshold` params.

**Files:**
- Modify: `nextflow_schema.json` — add to `pipeline_options` section

**Step 1: Add parameters to `pipeline_options.properties`**

Add after `netmhc3` (line 203):

```json
"run_rna_annotation": {
    "type": "boolean",
    "default": false,
    "description": "Enable RNA annotation of neoantigens using FORTE outputs.",
    "help_text": "When enabled, adds RNA expression (TPM) and RNA VAF columns to the neoantigen TSV output. Requires kallisto_abundance in samplesheet or rna_* columns in the input MAF.",
    "fa_icon": "fas fa-flask"
},
"run_fusion_neoantigens": {
    "type": "boolean",
    "default": false,
    "description": "Enable fusion neoantigen prediction from AGFusion output.",
    "help_text": "When enabled, processes AGFusion output from FORTE to predict fusion-derived neoantigens. Requires agfusion_dir in samplesheet.",
    "fa_icon": "fas fa-link"
},
"rna_tpm_threshold": {
    "type": "number",
    "default": 1.0,
    "description": "TPM threshold for classifying a gene as expressed.",
    "help_text": "Genes with TPM above this threshold (from Kallisto) are classified as expressed in the rna_expressed column.",
    "fa_icon": "fas fa-filter"
}
```

**Step 2: Commit**

```bash
git add nextflow_schema.json
git commit -m "feat: add RNA integration params to pipeline schema"
```

---

### Task 7: Update `nextflow.config` with Default Param Values

**Files:**
- Modify: `nextflow.config`

**Step 1: Find the params block and add defaults**

Add alongside existing params like `run_phylowgs`:

```groovy
run_rna_annotation      = false
run_fusion_neoantigens  = false
rna_tpm_threshold       = 1.0
```

**Step 2: Commit**

```bash
git add nextflow.config
git commit -m "feat: add default RNA param values to nextflow.config"
```

---

### Task 8: Update Main Workflow with RNA Annotation Channel

**Context:** This is the core pipeline integration. Add the optional RNA annotation path to `workflows/neoantigenpipeline.nf`.

**Files:**
- Modify: `workflows/neoantigenpipeline.nf`

**Step 1: Add module import**

At the top of the file (after line 19, the PHYLOWGS_STUB import):

```groovy
include { NEOANTIGENUTILS_RNAANNOTATE } from '../modules/msk/neoantigenutils/rnaannotate'
```

**Step 2: Parse optional samplesheet columns**

After the existing `ch_samplesheet` maps (around line 78), add:

```groovy
// Optional RNA channels
ch_maf_for_rna = ch_samplesheet.map { meta, maf, facets_hisens_cncf, hla_file ->
    [meta, maf]
}

ch_kallisto = ch_samplesheet.map { meta, maf, facets_hisens_cncf, hla_file, kallisto_abundance ->
    kallisto_abundance ? [meta, file(kallisto_abundance)] : [meta, []]
}
```

Note: The samplesheet map destructuring will need to be updated to handle the new optional columns. If the samplesheet has 4 columns (no RNA), the extra columns default to null. This depends on how nf-schema handles optional columns — check the existing `subworkflows/local/utils_nfcore_neoantigenpipeline_pipeline/main.nf` for the samplesheet parsing logic and update accordingly.

**Step 3: Add RNA annotation block after CONVERTANNOTJSON**

After line 131 (`NEOANTIGENUTILS_CONVERTANNOTJSON`), add:

```groovy
// Optional RNA annotation
if (params.run_rna_annotation) {
    NEOANTIGENUTILS_RNAANNOTATE(
        NEOANTIGENUTILS_CONVERTANNOTJSON.out.neoantigenTSV,
        ch_maf_for_rna,
        ch_kallisto,
        Channel.value(file(params.gtf))
    )
    ch_versions = ch_versions.mix(NEOANTIGENUTILS_RNAANNOTATE.out.versions)
}
```

**Step 4: Update emit block to include RNA outputs**

Add to the `emit:` section:

```groovy
rna_annotated_tsv = params.run_rna_annotation ?
    NEOANTIGENUTILS_RNAANNOTATE.out.annotated_tsv : Channel.empty()
rna_report        = params.run_rna_annotation ?
    NEOANTIGENUTILS_RNAANNOTATE.out.rna_report : Channel.empty()
```

**Step 5: Commit**

```bash
git add workflows/neoantigenpipeline.nf
git commit -m "feat: integrate RNA annotation into pipeline workflow"
```

---

### Task 9: Update Samplesheet Parsing for Optional Columns

**Context:** The samplesheet parsing in the utils subworkflow needs to handle the new optional columns gracefully.

**Files:**
- Modify: `subworkflows/local/utils_nfcore_neoantigenpipeline_pipeline/main.nf`

**Step 1: Examine current parsing logic**

Read the file and find where the samplesheet is parsed. The nf-schema plugin handles CSV parsing based on `schema_input.json`. Since the new columns are optional (not in `required`), rows without them will have null values.

**Step 2: Ensure the `ch_samplesheet` map handles 4 or 6 columns**

The existing map destructures `[meta, maf, facets_hisens_cncf, hla_file]`. With optional columns, nf-schema will pass null for missing values. Update the samplesheet map to pass through all columns:

```groovy
ch_samplesheet.map { row ->
    def meta = [id: row.sample]
    [meta, file(row.maf), file(row.facets_hisens_cncf), file(row.hla_file)]
}
```

The optional columns (`kallisto_abundance`, `agfusion_dir`) are accessed separately via their own channel maps that check for null.

**Step 3: Commit**

```bash
git add subworkflows/local/utils_nfcore_neoantigenpipeline_pipeline/main.nf
git commit -m "feat: update samplesheet parsing for optional RNA columns"
```

---

### Task 10: Add Fusion Neoantigen Path to Workflow

**Context:** The fusion path is more complex — it reuses NETMHCSTABANDPAN for MHC binding prediction on fusion peptides.

**Files:**
- Modify: `workflows/neoantigenpipeline.nf`

**Step 1: Add fusion module import**

```groovy
include { NEOANTIGENUTILS_FUSIONPREPARE } from '../modules/msk/neoantigenutils/fusionprepare'
```

**Step 2: Add fusion channel and processing block**

```groovy
// Optional fusion neoantigens
if (params.run_fusion_neoantigens) {
    ch_agfusion = ch_samplesheet.map { meta, maf, facets_hisens_cncf, hla_file, kallisto, agfusion_dir ->
        agfusion_dir ? [meta, file(agfusion_dir)] : null
    }.filter { it != null }

    NEOANTIGENUTILS_FUSIONPREPARE(ch_agfusion)
    ch_versions = ch_versions.mix(NEOANTIGENUTILS_FUSIONPREPARE.out.versions)

    // Fusion peptides go through existing MHC binding prediction
    ch_fusion_fasta_and_hla = NEOANTIGENUTILS_FUSIONPREPARE.out.mut_fasta
        .join(NEOANTIGENUTILS_FUSIONPREPARE.out.wt_fasta)
        .join(GENERATE_MUTATED_PEPTIDES.out.hla_string)

    // Note: Full fusion path through NETMHCSTABANDPAN -> NEOANTIGEN_EDITING
    // requires careful channel construction. This is the integration point
    // and may need adjustment based on how NETMHCSTABANDPAN handles
    // SV-only input vs mixed SNV+SV input.
}
```

**Step 3: Commit**

```bash
git add workflows/neoantigenpipeline.nf
git commit -m "feat: add fusion neoantigen processing path"
```

---

### Task 11: Add publishDir Configuration for New Outputs

**Files:**
- Modify: `conf/modules.config` (or equivalent config managing publishDir)

**Step 1: Add publishDir for RNA annotation outputs**

```groovy
withName: 'NEOANTIGENUTILS_RNAANNOTATE' {
    publishDir = [
        path: { "${params.outdir}/${meta.id}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename }
    ]
}

withName: 'NEOANTIGENUTILS_FUSIONPREPARE' {
    publishDir = [
        path: { "${params.outdir}/${meta.id}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename }
    ]
}
```

**Step 2: Commit**

```bash
git add conf/modules.config
git commit -m "feat: add publishDir config for RNA annotation outputs"
```

---

### Task 12: Create Test Data and Pipeline Test

**Context:** Create minimal test data for the RNA annotation path and add an nf-test.

**Files:**
- Create: `tests/data/test_abundance.tsv` (minimal Kallisto output)
- Create: `tests/data/test_rna_maf.maf` (MAF with rna_* columns)
- Modify: `tests/pipeline.nf.test` (add test case for RNA annotation)

**Step 1: Create test abundance file**

```tsv
target_id	length	eff_length	est_counts	tpm
ENST00000275493.7	3188	2939.42	1500	45.2
ENST00000269571.9	8560	8311.42	800	12.8
```

**Step 2: Create test MAF with RNA columns**

Copy existing `tests/data/dna.maf` and add `rna_t_alt_count`, `rna_t_ref_count`, `rna_t_variant_frequency` columns.

**Step 3: Add test config**

```groovy
// tests/rna_annotation.config
params {
    run_rna_annotation = true
    rna_tpm_threshold  = 1.0
}
```

**Step 4: Add nf-test case**

```groovy
test("Pipeline runs with RNA annotation") {
    when {
        params {
            input  = "${projectDir}/tests/data/rna_samplesheet.csv"
            outdir = "${outputDir}"
            run_rna_annotation = true
        }
    }
    then {
        assert workflow.success
        assert path("${outputDir}").list().any { it.contains("rna_annotated.tsv") }
    }
}
```

**Step 5: Commit**

```bash
git add tests/
git commit -m "test: add test data and nf-test for RNA annotation"
```

---

## Implementation Order Summary

| Task | What | Where | Dependencies |
|------|------|-------|-------------|
| 1 | `annotate_rna.py` + tests | Module repo | None |
| 2 | RNA annotation `main.nf` + meta | Module repo | Task 1 |
| 3 | `prepare_fusion_fasta.py` + tests | Module repo | None |
| 4 | Fusion prep `main.nf` + meta | Module repo | Task 3 |
| 5 | Update samplesheet schema | Pipeline repo | None |
| 6 | Update pipeline params schema | Pipeline repo | None |
| 7 | Update nextflow.config defaults | Pipeline repo | Task 6 |
| 8 | Add RNA annotation to workflow | Pipeline repo | Tasks 2, 5, 9 |
| 9 | Update samplesheet parsing | Pipeline repo | Task 5 |
| 10 | Add fusion path to workflow | Pipeline repo | Tasks 4, 8 |
| 11 | Add publishDir config | Pipeline repo | Tasks 8, 10 |
| 12 | Test data and pipeline tests | Pipeline repo | Tasks 8, 10 |

Tasks 1-4 (module repo) and Tasks 5-7 (pipeline repo) can be done in parallel.
Tasks 8-12 must be sequential after their dependencies.
