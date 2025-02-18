process VPT_SUM_SIGNALS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vpt:1.3' :
        'docker.io/vzgdocker/vpt:1.3' }"

    input:
    tuple val(meta), path(data)
    tuple val(meta), path(parquet)

    output:
    tuple val(meta), path("sum_signals.csv"), emit: signals
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def images = """${data}/images/mosaic_(?P<stain>[\\w|-]+)_z(?P<z>[0-9]+).tif"""
    def micron_to_mosaic = "${data}/images/micron_to_mosaic_pixel_transform.csv"


    """
    vpt \\
        --processes ${task.cpus} sum-signals \\
        --input-images="${images}" \\
        --input-boundaries ${parquet} \\
        --input-micron-to-mosaic ${micron_to_mosaic} \\
        --output-csv sum_signals.csv \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: 1.3
    END_VERSIONS
    """

}