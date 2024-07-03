> [![GitHub Actions CI Status](https://github.com/mskcc/neoantigenpipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/mskcc/neoantigenpipeline/actions/workflows/ci.yml) > [![GitHub Actions Linting Status](https://github.com/mskcc/neoantigenpipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/mskcc/neoantigenpipeline/actions/workflows/linting.yml) > [![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/mskcc/neoantigenpipeline)

## Introduction

**mskcc/neoantigenpipeline** is a bioinformatics pipeline that adapts Luksza et al.'s neoantigenEditing and fitness pipeline for usage by investigators in MSK. The pipeline curently supports working with TEMPO output mafs, Facets gene-level copy number calls, and Polysolver outputs. It outputs a json representation of the clonal structure of the tumor annotated with neoantigen burden, driver burden, and fitness of the clone. Also individual neoantigens are labeled with the quality of the neoantigen as described by Luksza et al.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Create phylogenetic trees using [PhyloWGS](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0602-8)
2. Use netMHCpan to calculate binding affinities [netMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
3. Use netMHCpanStab to calculate stability scores [netMHCpanStab](https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/)
4. Use Luksza et al.'s Neoantigen Quality and Fitness computations to evaluate peptides ([`Neoantigen Quality`](https://github.com/LukszaLab/NeoantigenEditing)))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,maf,facets_hisens_cncf,hla_file
tumor_normal,temp_test_somatic_unfiltered.maf,facets_hisens.cncf.txt,winners.hla.txt
tumor_normal2,temp_test_somatic_unfiltered.maf,facets_hisens.cncf.txt,winners.hla.txt
```
-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run mskcc/neoantigenpipeline \
   -profile prod,<docker/singularity> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

We thank the following people for their extensive assistance in the development of this pipeline:

- Nikhil ([@nikhil](https://github.com/nikhil))
- John ([@johnoooh](https://github.com/johnoooh))
- Alex ([@pintoa1-mskcc](https://github.com/pintoa1-mskcc))
- Martina ([@BradicM](https://github.com/BradicM))
- Allison ([@arichards2564](https://github.com/arichards2564))

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
