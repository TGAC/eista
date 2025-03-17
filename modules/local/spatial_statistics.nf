process SPATIAL_STATISTICS {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/villadsw/scanpy_squidpy:latest' :
        'docker.io/villadsw/scanpy_squidpy:latest' }"

    input:
    path h5ad_file
    // path samplesheet

    output:
    path "spatialstats"
    // path "spatialstats/*.h5ad",  emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    spatial_statistics.py \\
        --h5ad ${h5ad_file} \\
        --outdir spatialstats \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS        
    """



    // stub:
    // """
    // touch combined_matrix.h5ad
    // """
}
