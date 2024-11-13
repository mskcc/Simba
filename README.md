[![GitHub Actions CI Status](https://github.com/mskcc/simba/actions/workflows/ci.yml/badge.svg)](https://github.com/mskcc/simba/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/mskcc/simba/actions/workflows/linting.yml/badge.svg)](https://github.com/mskcc/simba/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/mskcc/simba)

## Introduction

**mskcc/simba** is a bioinformatics pipeline that runs the variant calling pipeline using other a collection of other nextflow workflows

1. Reallignment ([`Arrakis`](https://github.com/mskcc/arrakis))
2. Variant Calling ([`Odin`](https://github.com/mskcc/odin))
3. Structural Variants ([`Sif`](https://github.com/mskcc/sif))
4. Copy Number Variants ([`Loki`](https://github.com/mskcc/loki))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
pairId,tumorBam,normalBam,assay,normalType,bedFile
foo_sample,foo_tumor.rg.md.bam,foo_normal.rg.md.bam,IMPACT505,MATCHED,foo_tumor.foo_normal.fci.bed
bar_sample,bar_tumor.rg.md.bam,bar_normal.rg.md.bam,IMPACT505,MATCHED,bar_tumor.bar_normal.fci.bed
```

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run main.nf \
   -profile singularity,test_juno \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!NOTE]
> You must include `test_juno` as your profile for the workflow to run properly

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

mskcc/arrakis was originally written by Nikhil Kumar ([@nikhil](https://github.com/nikhil)).

We thank the following people for their extensive assistance in the development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

- Lisle E Mose, Charles M Perou, Joel S Parker, Improved indel detection in DNA and RNA via realignment with ABRA2, Bioinformatics, Volume 35, Issue 17, September 2019, Pages 2966–2973, https://doi.org/10.1093/bioinformatics/btz033
- “Picard Toolkit.” 2019. Broad Institute, GitHub Repository. https://broadinstitute.github.io/picard/; Broad Institute
- Van der Auwera, G. A., Carneiro, M. O., Hartl, C., Poplin, R., Del Angel, G., Levy-Moonshine, A., Jordan, T., Shakir, K., Roazen, D., Thibault, J., Banks, E., Garimella, K. V., Altshuler, D., Gabriel, S., & DePristo, M. A. (2013). From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Current protocols in bioinformatics, 43(1110), 11.10.1–11.10.33. https://doi.org/10.1002/0471250953.bi1110s43

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
