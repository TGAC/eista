process SPATIAL_TO_H5AD_XENIUM {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/villadsw/scanpy_squidpy_spatialdata:latest' :
        'docker.io/villadsw/scanpy_squidpy:latest' }"

    input:
    path outdir
    tuple val(meta), path(data)

    output:
    tuple val("raw"), path("*.h5ad"),  emit: h5ad
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tech = params.technology
    
    // run script
    """
    spatial_to_h5ad.py \\
        --tech ${tech} \\
        --datadir ${data} \\
        --outfile "${meta.id}_st_matrix.h5ad" \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """



    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/${meta.id}_matrix.h5ad
    touch versions.yml
    """
}
