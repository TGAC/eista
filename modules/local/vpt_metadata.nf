process VPT_METADATA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vpt:1.3' :
        'vzgdocker/vpt:1.3' }"

    input:
    tuple val(meta), path(parquet)
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("cell_metadata.csv"), emit: metadata
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    vpt \\
        --processes ${task.cpus} derive-entity-metadata \\
        --input-boundaries ${parquet} \\
        --input-entity-by-gene ${counts} \\
        --output-metadata cell_metadata.csv \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: 1.3
    END_VERSIONS
    """

}