/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_eista_pipeline'

include { VPT_QUANTIFICATION } from '../modules/local/vpt_quantification'
// include { VPT_SEGMENTATION } from '../modules/local/vpt_segmentation'
// include { VPT_PARTITION } from '../modules/local/vpt_partition'
// include { VPT_METADATA } from '../modules/local/vpt_metadata'
// include { VPT_SUM_SIGNALS } from '../modules/local/vpt_sum_signals'
// include { VPT_UPDATE_VZG } from '../modules/local/vpt_update_vzg'
// include { MTX_CONVERSION } from "../subworkflows/local/mtx_conversion"
include { SPATIAL_TO_H5AD } from '../modules/local/spatial_to_h5ad'
include { CONCAT_H5AD } from '../modules/local/concat_h5ad'
include { QC_CELL_FILTER } from "../modules/local/qc_cell_filter"
include { CLUSTERING_ANALYSIS } from "../modules/local/clustering_analysis"
include { SPATIAL_STATISTICS } from "../modules/local/spatial_statistics"
include { ANNOTATE_CELLS } from '../modules/local/annotate_cells'
include { TRAIN_CT_MODEL } from '../modules/local/train_ct_model'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EISTA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    samplesheet = file(params.input)
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_datacountsmeta = Channel.empty()
    ch_concat_h5ad_input = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     ch_samplesheet
    // )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())



    //===================================== Primary anaysis stage =====================================

    if (params.run_analyses.contains('primary')){

        if (params.technology == "vizgen") {
            if (!params.skip_analyses.contains('quantification')) {

                VPT_QUANTIFICATION (
                    ch_samplesheet,
                    ch_samplesheet.map{ meta, data -> file("${data[0]}/*.vzg") }
                )
                ch_versions = ch_versions.mix(VPT_QUANTIFICATION.out.versions)
                ch_datacountsmeta = VPT_QUANTIFICATION.out.datacountsmeta
                // ch_counts = VPT_QUANTIFICATION.out.counts
                // ch_metadata = VPT_QUANTIFICATION.out.metadata
            }
        }

        // Run mtx to h5ad conversion subworkflow
        // MTX_CONVERSION (
        //     ch_samplesheet,
        //     ch_counts,
        //     ch_metadata,
        //     samplesheet
        // )
        // ch_versions.mix(MTX_CONVERSION.out.ch_versions)

        // ch_samplesheet
        // .join(ch_counts, by: [0]).map{ meta, data, counts -> tuple(meta, data, counts) }
        // .join(ch_metadata, by: [0]).map{ meta, data, counts, metadata -> tuple(meta, data, counts, metadata) }
        // .set { ch_data_counts_meta }
        SPATIAL_TO_H5AD (
            Channel.value(file(params.outdir)),
            ch_datacountsmeta
        )
        // ch_versions = ch_versions.mix(SPATIAL_TO_H5AD.out.versions)

        ch_concat_h5ad_input = SPATIAL_TO_H5AD.out.h5ad.groupTuple()
        CONCAT_H5AD (
            ch_concat_h5ad_input,
            samplesheet
        )
        ch_versions = ch_versions.mix(CONCAT_H5AD.out.versions)  

    }

    //===================================== Secondary anaysis stage =====================================

    if (params.run_analyses.contains('secondary')){
    
        // MODULE: Run QC and cell filtering
        ch_h5ad = Channel.empty()
        if(params.run_analyses.contains('primary')){
            ch_h5ad = CONCAT_H5AD.out.h5ad
        }else if(params.h5ad){
            ch_h5ad = Channel.fromPath(params.h5ad)
        }else if(params.aligner){
            path = [
                'kallisto': "${params.outdir}/kallisto/mtx_conversions/combined_*_matrix.h5ad",
                'alevin': "${params.outdir}/alevin/mtx_conversions/combined_*_matrix.h5ad",
                'star': "${params.outdir}/star/mtx_conversions/combined_raw_matrix.h5ad"
            ].get(params.aligner)
            ch_h5ad = Channel.fromPath(path)
        }else{
            log.warn("For this analysis, please specify an h5ad file either by setting --aligner for the " +
            "h5ad file generated by the aligner or by setting --h5ad for an existing h5ad file.")
            return
        }

        if (!params.skip_analyses.contains('qccellfilter')) {
            QC_CELL_FILTER (
                ch_h5ad,
                Channel.fromPath(params.input)
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(QC_CELL_FILTER.out.versions)
            ch_h5ad = QC_CELL_FILTER.out.h5ad         
        }
        
        if (!params.skip_analyses.contains('clustering')) {
            CLUSTERING_ANALYSIS (
                ch_h5ad
            )
            ch_versions = ch_versions.mix(CLUSTERING_ANALYSIS.out.versions)
            ch_h5ad = CLUSTERING_ANALYSIS.out.h5ad
        } 

        if (!params.skip_analyses.contains('spatialstats')) {
            SPATIAL_STATISTICS (
                ch_h5ad
            )
            ch_versions = ch_versions.mix(SPATIAL_STATISTICS.out.versions)
        }          
        
    }


    //===================================== Tertiary anaysis stage =====================================

    if (!params.run_analyses.intersect(['tertiary', 'annotation', 'dea']).is()){
    
        // Get input h5ad file
        ch_h5ad = Channel.empty()
        if(params.run_analyses.contains('secondary')){
            if (!params.skip_analyses.contains('clustering')) {
                ch_h5ad = CLUSTERING_ANALYSIS.out.h5ad
            }else {
                ch_h5ad = QC_CELL_FILTER.out.h5ad
            }
        }else if(params.h5ad){
            ch_h5ad = Channel.fromPath(params.h5ad)
        }else{
            path1 = "${params.outdir}/clustering/adata_clustering.h5ad"
            path2 = "${params.outdir}/qc_cell_filter/adata_filtered_normalized.h5ad"
            path3 = "${params.outdir}/annotation/adata_annotation.h5ad"
            if(params.run_analyses.contains('dea') && (new File(path3).exists())){
                ch_h5ad = Channel.fromPath(path3)           
            }else if(new File(path1).exists()){
                ch_h5ad = Channel.fromPath(path1)
            }else if(new File(path2).exists()){
                ch_h5ad = Channel.fromPath(path2)
            }
        }
        ch_h5ad.ifEmpty {
            log.warn("For this analysis, h5ad file can be found in secondary analysis, please specify an h5ad file by setting --h5ad.")
            return            
        }

        if (params.run_analyses.any{it=='tertiary' || it=='annotation'} and !params.skip_analyses.contains('annotation')) {
            // ch_ctmodel = params.ctmodel? Channel.fromPath(params.ctmodel) : []
            ANNOTATE_CELLS (
                ch_h5ad,
                // ch_ctmodel
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(ANNOTATE_CELLS.out.versions)
            ch_h5ad = ANNOTATE_CELLS.out.h5ad      
        }
        
        // if (params.run_analyses.any{it=='tertiary' || it=='dea'} and !params.skip_analyses.contains('dea')) {
        //     RANK_GENES (
        //         ch_h5ad,
        //     )
        //     ch_versions = ch_versions.mix(RANK_GENES.out.versions)
        // }
   
        
    }

    // train cell-type models for CellTypist
    if (params.run_analyses.contains('ctmodel') and !params.skip_analyses.contains('ctmodel')) {
            ch_h5ad = Channel.fromPath(params.h5ad)
            TRAIN_CT_MODEL (
                ch_h5ad,
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(TRAIN_CT_MODEL.out.versions)
            // ch_h5ad = ANNOTATE_CELLS.out.h5ad      
    }





    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

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

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
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
        ch_multiqc_logo.toList()
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
