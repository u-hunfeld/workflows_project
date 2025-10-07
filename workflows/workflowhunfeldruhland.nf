/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { CUTADAPT               } from '../modules/nf-core/cutadapt/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_workflowhunfeldruhland_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WORKFLOWHUNFELDRUHLAND {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Handle genome parameters - assign from iGenomes if genome is specified
    //
    def fasta_file = null
    def gtf_file = null
    def star_index_dir = null

    if (params.genome && params.genomes && params.genomes.containsKey(params.genome)) {
        fasta_file = params.genomes[params.genome].fasta
        gtf_file = params.genomes[params.genome].gtf
        star_index_dir = params.genomes[params.genome].star
        log.info "Using iGenomes for genome: ${params.genome}"
    } else {
        fasta_file = params.fasta
        gtf_file = params.gtf
        star_index_dir = params.star_index
        log.info "Using manual genome parameters"
    }

    // Debug output
    log.info "Using FASTA: ${fasta_file}"
    log.info "Using GTF: ${gtf_file}"
    log.info "Using STAR index: ${star_index_dir}"

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run CUTADAPT
    //
    CUTADAPT(
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
    
    // Use trimmed reads for downstream analysis
    ch_trimmed_reads = CUTADAPT.out.reads

    ch_trimmed_reads.view() //for debugging, can be deleted again later

    //
    // Prepare STAR index and reference files
    //
    if (star_index_dir) {
        // Use existing STAR index
        ch_star_index = Channel.fromPath(star_index_dir, checkIfExists: true)
        ch_gtf = Channel.fromPath(gtf_file, checkIfExists: true)
    } else if (fasta_file && gtf_file) {
        // Generate STAR index from FASTA and GTF
        
        ch_fasta= Channel.fromPath(fasta_file).map{fasta -> [[:], fasta]}
        ch_gtf_for_index=Channel.fromPath(gtf_file).map{gtf -> [[:], gtf]}

        STAR_GENOMEGENERATE (
            ch_fasta, ch_gtf_for_index
        )
        ch_star_index = STAR_GENOMEGENERATE.out.index
        ch_gtf = Channel.fromPath(gtf_file, checkIfExists: true)
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    } else {
        error "Either provide --star_index, or both --fasta and --gtf, or use --genome to auto-select reference files"
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_workflowhunfeldruhland_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: Align reads with STAR
    //
    STAR_ALIGN (
        ch_trimmed_reads,
        ch_star_index,
        ch_gtf,
        false,                    // star_ignore_sjdbgtf (boolean, not channel)
        params.seq_platform ?: '',  // string value, not channel
        params.seq_center ?: ''     // string value, not channel
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // Collect MultiQC files
    //
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]})

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
