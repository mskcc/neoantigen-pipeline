process PHYLOWGS_WRITERESULTSGO {
    tag "$meta.id"
    label 'process_medium'
    // Reuses the existing phylowgs container: write_results_go.py (bin/) only
    // replaces the pickle-based tree loader, it still calls the same
    // pwgsresults.ResultMunger / JsonWriter classes this image already ships.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/mskcc-omics-workflows/phylowgs:v1.5-msk':
        'ghcr.io/mskcc-omics-workflows/phylowgs:v1.5-msk' }"

    input:
    tuple val(meta), path(trees), path(mutlist)

    output:
    tuple val(meta), path("*.summ.json.gz")     , emit: summ
    tuple val(meta), path("*.muts.json.gz")     , emit: muts
    tuple val(meta), path("*.mutass.zip")       , emit: mutass
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    write_results_go.py \\
        ${args} \\
        --include-ssm-names \\
        ${prefix} \\
        ${trees} \\
        ${mutlist} \\
        ${prefix}.summ.json.gz \\
        ${prefix}.muts.json.gz \\
        ${prefix}.mutass.zip
    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    phylowgs: \$PHYLOWGS_TAG
	END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.summ.json.gz
    touch ${prefix}.muts.json.gz
    touch ${prefix}.mutass.zip

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    phylowgs: \$PHYLOWGS_TAG
	END_VERSIONS
    """
}
