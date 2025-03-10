process VPT_METADATA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vpt:1.3' :
        'docker.io/vzgdocker/vpt:1.3' }"

    input:
    tuple val(meta), path(parquet), path(counts)
    // tuple val(meta), path(counts)

    output:
    tuple val(meta), path("cell_metadata.csv"), emit: metadata
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export HOME=$PWD
    export CELLPOSE_LOCAL_MODELS_PATH=/ei/projects/0/05407428-a659-41d6-a7bd-4567bf45e494/data/vizgen/models
    vpt \\
        --verbose derive-entity-metadata \\
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