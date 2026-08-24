process NEOANTIGENUTILS_NEOANTIGENINPUT {
    tag "$meta.id"
    label 'process_medium'
    // TEMPORARY: the plain :1.6.1 tag's multi-arch manifest is stale (missing the
    // unzip fix) because create-ghcr-manifest was skipped upstream when an unrelated
    // gbcms push failed in the same CI matrix run. The per-arch :1.6.1-amd64 tag was
    // built and pushed correctly and already has unzip. Revert to :1.6.1 once the
    // manifest is re-stitched (mskcc-omics-workflows/containers).
    container "ghcr.io/mskcc-omics-workflows/neoantigen-utils-base:1.6.1-amd64"

    input:
    tuple val(meta),  path(inputMaf),      path(inputBedpe, arity: '0..*'),    path(hlaFile)
    tuple val(meta2), path(phyloWGSsumm),  path(phyloWGSmut),   path(phyloWGSfolder)
    tuple val(meta3), path(mutNetMHCpan),  path(wtNetMHCpan)
    tuple path(gtf),  path(cdna)

    output:
    tuple val(meta), path("*_input.json"),                                                  emit: json
    path "versions.yml",                                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def id = task.ext.prefix ?: "${meta.id}"
    def patientid = task.ext.cohort ?: "${meta.id}_patient"
    def cohort = task.ext.cohort ?: "${meta.id}_cohort"
    def bedpe = inputBedpe ? "--bedpe_file ${inputBedpe}": ""

    """
        # TEMPORARY: pyensembl (via generate_input.py's ensembl_load) defaults its
        # index cache to \$HOME, which is read-only on iris compute nodes. Redirect
        # it into the task's own writable work dir. Revert if generate_input.py
        # grows a --cache-dir flag or pyensembl gains a container-friendly default.
        export PYENSEMBL_CACHE_DIR="\$PWD/.pyensembl_cache"

        tree_folder_name=\$(basename -s .zip "${phyloWGSfolder}")
        mkdir \$tree_folder_name
        unzip ${phyloWGSfolder} -d \$tree_folder_name
        gzip -d -c ${phyloWGSsumm} > ${id}.summ.json
        gzip -d -c ${phyloWGSmut} > ${id}.mut.json



        generate_input.py --maf_file ${inputMaf} \
        ${bedpe} \
        --summary_file ${id}.summ.json \
        --mutation_file ${id}.mut.json \
        --tree_directory \$tree_folder_name \
        --id ${id} --patient_id ${patientid} \
        --cohort ${cohort} --HLA_genes ${hlaFile} \
        --netMHCpan_MUT_input ${mutNetMHCpan} \
        --netMHCpan_WT_input ${wtNetMHCpan} \
        --gtf-file ${gtf} \
        --cdna-file ${cdna} \
        ${args}

        cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    neoantigeninput: \$(echo \$(generate_input.py -v))
	END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def id = task.ext.prefix ?: "${meta.id}"
    def patientid =task.ext.cohort ?: "${meta.id}_patient"
    def cohort =task.ext.cohort ?: "${meta.id}_cohort"
    """

        touch ${patientid}_${id}_input.json

        cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    neoantigeninput: \$(echo \$(generate_input.py -v))
	END_VERSIONS
    """
}
