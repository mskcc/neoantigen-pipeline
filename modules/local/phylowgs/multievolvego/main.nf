process PHYLOWGS_MULTIEVOLVEGO {
    tag "$meta.id"
    label 'process_high'
    // Image built from the go-port branch Dockerfile at mskcc/phylowgs (PhyloWGS_refactor/Dockerfile).
    // Not yet published to ghcr — temporarily published to orgeraj/phylowgs-go on Docker Hub
    // (multi-arch: linux/amd64 + linux/arm64) for cohort testing until the ghcr image lands.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://orgeraj/phylowgs-go:go-port':
        'orgeraj/phylowgs-go:go-port' }"

    input:
    tuple val(meta), path(cnv_data), path(ssm_data)

    output:
    tuple val(meta), path("chains/trees.zip"), path("chains/mutlist.json"), emit: trees
    tuple val(meta), path("chains/best_tree.json"), emit: best_tree
    tuple val(meta), path("chains/mutlist.json")  , emit: mutlist
    tuple val(meta), path("chains/summary.json")  , emit: summary
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    phylowgs-go \\
        ${args} \\
        -O chains \\
        -D ${prefix} \\
        ${ssm_data} \\
        ${cnv_data}

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    phylowgs_go: \$PHYLOWGS_GO_TAG
	END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir chains
    touch chains/trees.zip chains/best_tree.json chains/mutlist.json chains/summary.json

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    phylowgs_go: \$PHYLOWGS_GO_TAG
	END_VERSIONS
    """
}
