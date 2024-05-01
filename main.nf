#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mskcc/simba
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/mskcc/simba
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ARRAKIS  } from './workflows/arrakis/workflows/arrakis.nf'
include { CNV } from './workflows/loki/subworkflows/local/cnv'
include { SV } from './workflows/sif/subworkflows/local/sv'
include { CALL_VARIANTS } from './workflows/odin/subworkflows/local/variant-calling/main'
include { FIND_COVERED_INTERVALS } from './workflows/odin/subworkflows/local/find_covered_intervals'
include { MAF_PROCESSING } from './workflows/odin/subworkflows/local/maf-processing/main'
include { MAF_FILTER_WORKFLOW } from './workflows/odin/subworkflows/local/maf-filter/main'
include { MAF_ANNOTATE } from './workflows/odin/subworkflows/local/maf-annotate/main'
include { TMB_WORKFLOW } from './workflows/odin/subworkflows/local/tmb/main'
include { SAMTOOLS_HEADER_VIEW as normal_header; SAMTOOLS_HEADER_VIEW as tumor_header} from './workflows/odin/modules/local/get_bam_header'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_simba_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_simba_pipeline'
include { softwareVersionsToYAML } from './subworkflows/nf-core/utils_nfcore_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_simba_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow MSKCC_SIMBA {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SIMBA (
        samplesheet
    )

    emit:
    multiqc_report = SIMBA.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    // Import resources

    ch_fasta_ref = Channel.value([ "reference_genome", file(params.fasta) ])
    ref_index_list = []
    for(single_genome_ref in params.fasta_index){
        ref_index_list.add(file(single_genome_ref))
    }
    ch_fasta_fai_ref = Channel.value([ "reference_genome_index",ref_index_list])
    ch_delly_exclude = Channel.value([ "delly_exclude", file(params.delly_exclude) ])
    ch_exac_filter = Channel.value(["exac_filter", file(params.exac_filter)])
    ch_exac_filter_index = Channel.value(["exac_filter_index", file(params.exac_filter_index)])
    ch_bedfile = Channel.value([ file(params.bed_file) ])
    ch_dbsnp = Channel.value([ "dbsnp", file(params.dbsnp) ])
    ch_cosmic = Channel.value([ "cosmic", file(params.cosmic) ])
    ch_hotspot = Channel.value([ "hotspot", file(params.hotspot) ])
    ch_impact_gene_list = Channel.value([ "impact_gene_list", file(params.impact_gene_list) ])
    ch_curated_bams = get_curated_bams()

    //
    // WORKFLOW: Run main workflow
    //

    ch_versions = Channel.empty()



    ARRAKIS (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    normal_header (
        ARRAKIS.out.normal_bam
            .map{
                new Tuple(it[0],it[1])
            }
    )

    tumor_header (
        ARRAKIS.out.tumor_bam
            .map{
                new Tuple(it[0],it[1])
            }
    )

    ch_versions = ch_versions.mix(ARRAKIS.out.versions)

    // Run Loki workflow



    input_bam_channel = create_bam_channel(ARRAKIS.out.normal_bam, ARRAKIS.out.tumor_bam)
    input_bam_channel_with_bed = create_bam_channel_with_sample_name(ARRAKIS.out.normal_bam, ARRAKIS.out.tumor_bam,normal_header.out.sample_name,tumor_header.out.sample_name)

    ch_normal = input_bam_channel_with_bed
        .map{
            new Tuple(it[0],it[1][1], it[2][1])
        }
    ch_tumor = input_bam_channel_with_bed
        .map{
            new Tuple(it[0],it[1][0], it[2][0])
        }

    CNV (
        input_bam_channel
    )

    ch_versions = ch_versions.mix(CNV.out.versions)

    SV (
        ch_normal,
        ch_tumor,
        ch_fasta_ref,
        ch_fasta_fai_ref,
        ch_delly_exclude,
        params.delly_type,
        ch_exac_filter,
        ch_exac_filter_index
    )

    ch_versions = ch_versions.mix(SV.out.versions)

    intervals = params.intervals

    input_files = input_bam_channel_with_bed
        .branch {
            no_bed: it[3] == null
            bed_provided: it[3] != null
        }

    fci_input = input_files.no_bed
        .map{
            new Tuple(it[0],it[1],it[2])
        }

    FIND_COVERED_INTERVALS (
        fci_input,
        ch_fasta_ref,
        ch_fasta_fai_ref,
        intervals
    )

    ch_versions = ch_versions.mix(FIND_COVERED_INTERVALS.out.versions)

    generated_beds = join_bams_with_bed(fci_input, FIND_COVERED_INTERVALS.out.bed_file)

    variant_input = input_files.bed_provided.mix(generated_beds)

    CALL_VARIANTS (
        variant_input,
        ch_fasta_ref,
        ch_fasta_fai_ref,
        ch_dbsnp,
        ch_cosmic,
        ch_hotspot
    )

    ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)

    MAF_PROCESSING (
        CALL_VARIANTS.out.annotate_vcf,
        ch_fasta_ref,
        ch_fasta_fai_ref,
        ch_exac_filter,
        ch_exac_filter_index,
        input_bam_channel,
        ch_curated_bams
    )

    ch_versions = ch_versions.mix(MAF_PROCESSING.out.versions)

    MAF_FILTER_WORKFLOW (
        MAF_PROCESSING.out.maf
    )

    ch_versions = ch_versions.mix(MAF_FILTER_WORKFLOW.out.versions)

    MAF_ANNOTATE (
        MAF_FILTER_WORKFLOW.out.analysis_maf,
        ch_impact_gene_list
    )

    ch_versions = ch_versions.mix(MAF_ANNOTATE.out.versions)

    TMB_WORKFLOW (
        MAF_FILTER_WORKFLOW.out.data_mutations_extended_file
    )

    ch_versions = ch_versions.mix(TMB_WORKFLOW.out.versions)

    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    //params.input = "${params.outdir}/realignment_bams_samplesheet.csv"

    //LOKI {
    //    ARRAKIS.out.samplesheet
    //}



    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url
    )
}

def create_bam_channel(normal, tumor) {
    tumor_channel = tumor
        .map{
            new Tuple(it[0].id,it)
            }
    normal_channel = normal
        .map{
            new Tuple(it[0].id,it)
            }
    mergedWithKey = tumor_channel
        .join(normal_channel)
    merged = mergedWithKey
        .map{
            new Tuple(it[1][0], it[2][1],it[2][2], it[1][1],it[1][2],[],[])
        }
    return merged
}

def join_bams_with_bed(bams,bed) {
        bam_channel = bams
            .map{
                new Tuple(it[0].id,it)
                }
        bed_channel = bed
            .map{
                new Tuple(it[0].id,it)
                }
        mergedWithKey = bam_channel
            .join(bed_channel)
        merged = mergedWithKey
            .map{
                new Tuple(it[1][0],it[1][1],it[1][2],it[2][1])
            }
        return merged

}

// Gets a list of curated bams
def get_curated_bams() {

    bams = files("${params.curated_bam_path}/*.bam")
    bam_indexes = files("${params.curated_bam_path}/*.bai")

    if(bams.size == 0)
    {
        exit 1, "ERROR: Please check the curated bams folder contains bams: \n${params.curated_bam_path}"
    }


    if(bam_indexes.size == 0)
    {
        exit 1, "ERROR: Please check the curated bams folder contains bai index files: \n${params.curated_bam_path}"

    }

    curated_bams = Channel.value(["curated_bams", bams, bam_indexes])
    return curated_bams


}


def create_bam_channel_with_sample_name(normal, tumor, normal_sample_name, tumor_sample_name) {
    tumor_channel = tumor
        .map{
            new Tuple(it[0].id,it)
            }
    tumor_name_channel = tumor_sample_name
        .map{
            new Tuple(it[0].id,it)
            }
    normal_channel = normal
        .map{
            new Tuple(it[0].id,it)
            }
    normal_name_channel = normal_sample_name
        .map{
            new Tuple(it[0].id,it)
            }
    mergedWithKey = tumor_channel
        .join(normal_channel)
        .join(tumor_name_channel)
        .join(normal_name_channel)
    merged = mergedWithKey
        .map{
            meta = it[1][0]
            normalSample = it[4][1]
            tumorSample = it[3][1]
            meta.tumorSampleName = tumorSample.trim()
            meta.normalSampleName = normalSample.trim()
            new Tuple(it[1][0], [it[1][1],it[2][1]], [it[1][2],it[2][2]],null)
        }

    return merged
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
