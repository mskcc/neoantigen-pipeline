process NEOANTIGENEDITING_COMPUTEFITNESS {
    tag "$meta.id"
    label 'process_medium'

    // TEMPORARY: the plain :1.1 tag's multi-arch manifest is stale (missing the
    // PATH/exec-bit fix for compute_fitness.py) for the same reason as
    // neoantigen-utils-base:1.6.1 -- create-ghcr-manifest was skipped upstream when
    // an unrelated gbcms push failed in the same CI matrix run. The per-arch
    // :1.1-amd64/-arm64 tags built and pushed correctly. Revert to :1.1 once the
    // manifest is re-stitched (mskcc-omics-workflows/containers).
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/mskcc-omics-workflows/neoantigen-editing:1.1-amd64':
        'ghcr.io/mskcc-omics-workflows/neoantigen-editing:1.1-amd64' }"

    input:
    tuple val(meta),  path(patient_data), path(alignment_file)

    output:
    tuple val(meta), path("*_annotated.json")               , emit: annotated_output
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    compute_fitness.py \\
        --alignment ${alignment_file} \\
        --sample_file ${patient_data} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    neoantigenEditing: \$NEOANTIGEN_EDITING_TAG
	END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch patient_data_annotated.json

    cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    neoantigenEditing: \$NEOANTIGEN_EDITING_TAG
	END_VERSIONS
    """
}
