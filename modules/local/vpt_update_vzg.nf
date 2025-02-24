process VPT_UPDATE_VZG {
    tag "$meta.id"
    label 'process_process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vpt:1.3' :
        'docker.io/vzgdocker/vpt:1.3' }"

    input:
    tuple val(meta), path(data)
    tuple val(meta), path(parquet)
    tuple val(meta), path(counts)
    tuple val(meta), path(metadata)

    output:
    tuple val(meta), path($vzg_filename), emit: vzg
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vzg = "${data}/*.vzg"
    def vzg_filename = 'update_' + new File(vzg).getName()


    """
    if [ -f ${vzg} ]; then
        vpt \\
            --processes ${task.cpus} update-vzg \\
            --input-vzg ${vzg} \\
            --input-boundaries ${parquet} \\
            --input-entity-by-gene ${counts} \\
            --input-metadata ${metadata} \\
            --output-csv ${vzg_filename} \\
            $args \\
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: 1.3
    END_VERSIONS
    """

}