/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { SPATIAL_TO_H5AD   }             from '../../modules/local/spatial_to_h5ad.nf'
include { CONCAT_H5AD   }             from '../../modules/local/concat_h5ad.nf'
// include { MTX_TO_SEURAT }             from '../../modules/local/mtx_to_seurat.nf'

workflow MTX_CONVERSION {

    take:
    data
    counts
    metadata
    samplesheet

    main:
        ch_versions = Channel.empty()

        //
        // Convert matrix to h5ad
        //
        SPATIAL_TO_H5AD (
            // data,
            counts,
            metadata
        )

        //
        // Concat sample-specific h5ad in one
        //
        ch_concat_h5ad_input = SPATIAL_TO_H5AD.out.h5ad.collect() // gather all sample-specific files / per type
        CONCAT_H5AD (
            ch_concat_h5ad_input,
            samplesheet
        )

        //
        // Convert matrix do seurat
        //
        // MTX_TO_SEURAT (
        //     mtx_matrices
        // )

        //TODO CONCAT h5ad and MTX to h5ad should also have versions.yaml output
        ch_versions = ch_versions.mix(SPATIAL_TO_H5AD.out.versions)
        // ch_versions = ch_versions.mix(SPATIAL_TO_H5AD.out.versions, MTX_TO_SEURAT.out.versions)

    emit:
    ch_versions
    h5ad = CONCAT_H5AD.out.h5ad
    // counts = MTX_TO_H5AD.out.counts  was this ever used?

}
