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
    // def images = "${data}/images/mosaic_*_z*.tif"
    def transformation = "${data}/images/micron_to_mosaic_pixel_transform.csv"
    // if(!new File(transformation).exists()) {transformation = "${data}/micron_to_mosaic_pixel_transform.csv"}

    """
    export HOME=$PWD
    export CELLPOSE_LOCAL_MODELS_PATH=/ei/projects/0/05407428-a659-41d6-a7bd-4567bf45e494/data/vizgen/models
    vpt \\
        --verbose --processes ${task.cpus} run-segmentation \\
        --segmentation-algorithm ${algorithm} \\
        --input-images="${images}" \\
        --input-micron-to-mosaic ${transformation} \\
        --output-path "./" \\
        --tile-size ${params.tile_size} \\
        --tile-overlap ${params.tile_overlap} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: 1.3
    END_VERSIONS
    """

}