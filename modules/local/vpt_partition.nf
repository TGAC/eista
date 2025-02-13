process VPT_PARTITION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vpt:1.3' :
        'vzgdocker/vpt:1.3' }"

    input:
    tuple val(meta), path(parquet)

    output:
    tuple val(meta), path("cell_by_gene.csv"), emit: counts
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def transcripts = new File("${data}/detected_transcripts.csv").text


    """
    vpt \\
        --processes ${task.cpus} partition-transcripts \\
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