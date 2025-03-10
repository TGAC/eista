process VPT_PARTITION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vpt:1.3' :
        'docker.io/vzgdocker/vpt:1.3' }"

    input:
    tuple val(meta), path(data), path(parquet)
    // tuple val(meta), path(parquet)

    output:
    tuple val(meta), path("cell_by_gene.csv"), emit: counts
    tuple val(meta), path("detected_transcripts.csv"), emit: transcripts
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def transcripts = "${data}/detected_transcripts.csv"


    """
    export HOME=$PWD
    export CELLPOSE_LOCAL_MODELS_PATH=/ei/projects/0/05407428-a659-41d6-a7bd-4567bf45e494/data/vizgen/models
    vpt \\
        --verbose partition-transcripts \\
        --input-boundaries ${parquet} \\
        --input-transcripts ${transcripts} \\
        --output-entity-by-gene cell_by_gene.csv \\
        --output-transcripts detected_transcripts.csv \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: 1.3
    END_VERSIONS
    """

}