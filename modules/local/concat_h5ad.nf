process CONCAT_H5AD {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/villadsw/scanpy_squidpy:latest' :
        'docker.io/villadsw/scanpy_squidpy:latest' }"

    input:
    path(h5ad)
    path samplesheet

    output:
    path "*.h5ad", emit: h5ad
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    concat_h5ad.py \\
        --input $samplesheet \\
        --out combined_st_matrix.h5ad \\
        --suffix "_st_matrix.h5ad"
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    // --suffix "_matrix.h5ad" : change to remove 'input_type' so that obs can map samplesheet

    stub:
    """
    touch combined_matrix.h5ad
    """
}
