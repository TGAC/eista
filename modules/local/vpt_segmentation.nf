process VPT_SEGMENTATION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vpt:1.3' :
        'docker.io/vzgdocker/vpt:1.3' }"

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("cellpose_micron_space.parquet"), emit: parquet
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def algorithm = (params.segmentation_algorithm)?: new File("${workflow.projectDir}/assets/algorithms/cellpose_default_1_ZLevel.json").text
    def algorithm = (params.segmentation_algorithm)?: "${workflow.projectDir}/assets/algorithms/cellpose_default_1_ZLevel.json"
    // def images = new File("${data}/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif").text
    def images = """${data}/images/mosaic_(?P<stain>[\\w|-]+)_z(?P<z>[0-9]+).tif"""
    def micron_to_mosaic = "${data}/images/micron_to_mosaic_pixel_transform.csv"
    // if(!new File(micron_to_mosaic).exists()) {micron_to_mosaic = "${data}/micron_to_mosaic_pixel_transform.csv"}

    """
    vpt \\
        --processes ${task.cpus} run-segmentation \\
        --segmentation-algorithm ${algorithm} \\
        --input-images="${images}" \\
        --input-micron-to-mosaic ${micron_to_mosaic} \\
        --output-path "./" \\
        --tile-size ${params.tile_size} \\
        --tile-overlap ${params.tile_overlap} \\
        --overwrite \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: 1.3
    END_VERSIONS
    """

}