process PHYLOWGS_STUB {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(maf), path(facets_genelevel)

    output:
    tuple val(meta), path("*.summ.json.gz"),  emit: summ
    tuple val(meta), path("*.muts.json.gz"),  emit: muts
    tuple val(meta), path("*.mutass.zip"),    emit: mutass
    path "versions.yml",                      emit: versions

    script:
    def id = meta.id
    """
    echo '{"trees": {}, "params": {}, "dataset_name": "", "tree_densities": {}}' | gzip > ${id}.summ.json.gz
    echo '{"ssms": {}}' | gzip > ${id}.muts.json.gz
    mkdir empty_tree_dir && zip -r ${id}.mutass.zip empty_tree_dir
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phylowgs_stub: "1.0"
    END_VERSIONS
    """
}
