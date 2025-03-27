process SPATIAL_TO_H5AD {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/villadsw/scanpy_squidpy:latest' :
        'docker.io/villadsw/scanpy_squidpy:latest' }"

    input:
    path outdir
    tuple val(meta), path(counts), path(metadata)
    // tuple val(meta), path(counts)
    // tuple val(meta), path(metadata)

    output:
    tuple val("raw"), path("*.h5ad"),  emit: h5ad
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tech = params.technology
    // def transformation = "${data}/images/micron_to_mosaic_pixel_transform.csv"
    def datadir = "${outdir}/${params.technology}/${meta.id}"
    def transformation 
    if (params.technology =='vizgen') {
        
        // transformation = "${data}/images/micron_to_mosaic_pixel_transform.csv"
        transformation = "micron_to_mosaic_pixel_transform.csv"

    } else if (params.aligner == 'kallisto') {

        kb_pattern   = (input_type == 'raw') ? 'un' : ''
        mtx_dir      = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : "counts_${kb_pattern}filtered"
        if ((input_type == 'custom_emptydrops_filter') && (params.kb_workflow != 'standard')) { mtx_dir = 'emptydrops_filtered/\${input_type}' } // dir has subdirs for non-standard workflows
        mtx_matrix   = "${mtx_dir}/*.mtx"
        barcodes_tsv = "${mtx_dir}/*.barcodes.txt"
        features_tsv = "${mtx_dir}/*.genes.names.txt"

        kb_non_standard_files = ""
        if (params.kb_workflow == "lamanno") {
            kb_non_standard_files = "spliced unspliced"
            matrix       = "${mtx_dir}/\${input_type}.mtx"
            barcodes_tsv = "${mtx_dir}/\${input_type}.barcodes.txt"
            features_tsv = "${mtx_dir}/\${input_type}.genes.txt"
        }
        if (params.kb_workflow == "nac") {
            kb_non_standard_files = "nascent ambiguous mature"
            matrix       = "${mtx_dir}/*\${input_type}.mtx"
            features_tsv = "${mtx_dir}/*.genes.txt"
        }

    } else if (params.aligner == 'alevin') {

        // alevin does not have filtered/unfiltered results
        mtx_dir      = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : '*_alevin_results/af_quant/alevin'
        mtx_matrix   = "${mtx_dir}/quants_mat.mtx"
        barcodes_tsv = "${mtx_dir}/quants_mat_rows.txt"
        features_tsv = "${mtx_dir}/quants_mat_cols.txt"

    } else if (params.aligner == 'star') {

        mtx_dir      = (input_type == 'custom_emptydrops_filter') ? 'emptydrops_filtered' : "${input_type}"
        suffix       = (input_type == 'custom_emptydrops_filter') ? '' : '.gz'
        mtx_matrix   = "${mtx_dir}/matrix.mtx${suffix}"
        barcodes_tsv = "${mtx_dir}/barcodes.tsv${suffix}"
        features_tsv = "${mtx_dir}/features.tsv${suffix}"

    }

    //
    // run script
    //
    if (params.technology =='vizgen')
    """
    # convert file types
    spatial_to_h5ad.py \\
        --tech ${tech} \\
        --datadir ${datadir} \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --transformation micron_to_mosaic_pixel_transform.csv \\
        --outfile "${meta.id}_st_matrix.h5ad" \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    else if (params.aligner == 'kallisto' && params.kb_workflow != 'standard')
    """
    # convert file types
    for input_type in ${kb_non_standard_files} ; do
        spatial_to_h5ad.py \\
            --aligner ${params.aligner} \\
            --sample ${meta.id} \\
            --input ${matrix} \\
            --barcode ${barcodes_tsv} \\
            --feature ${features_tsv} \\
            --txp2gene ${txp2gene} \\
            --star_index ${star_index} \\
            --out ${meta.id}/${meta.id}_\${input_type}_matrix.h5ad ;
    done
    """

    else
    """
    # convert file types
    spatial_to_h5ad.py \\
        --task_process ${task.process} \\
        --aligner ${params.aligner} \\
        --sample ${meta.id} \\
        --input $mtx_matrix \\
        --barcode $barcodes_tsv \\
        --feature $features_tsv \\
        --txp2gene ${txp2gene} \\
        --star_index ${star_index} \\
        --out ${meta.id}/${meta.id}_${input_type}_matrix.h5ad
    """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/${meta.id}_matrix.h5ad
    touch versions.yml
    """
}
