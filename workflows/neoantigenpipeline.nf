/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_neoantigenpipeline_pipeline'
include { PHYLOWGS_CREATEINPUT } from '../modules/msk/phylowgs/createinput/main'
include { PHYLOWGS_MULTIEVOLVE } from '../modules/msk/phylowgs/multievolve/main'
include { PHYLOWGS_PARSECNVS } from '../modules/msk/phylowgs/parsecnvs/main'
include { PHYLOWGS_WRITERESULTS } from '../modules/msk/phylowgs/writeresults/main'
include { PHYLOWGS } from '../subworkflows/msk/phylowgs'
include { NEOANTIGENINPUT } from '../modules/msk/neoantigeninput/main'
include { NETMHCPAN } from '../modules/msk/netmhcpan/main'
include { NEOANTIGENEDITING_ALIGNTOIEDB } from '../modules/msk/neoantigenediting/aligntoiedb'
include { NEOANTIGENEDITING_COMPUTEFITNESS } from '../modules/msk/neoantigenediting/computefitness'
include { NEOANTIGEN_EDITING } from '../subworkflows/msk/neoantigen_editing'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NEOANTIGENPIPELINE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input It should have maf, polysolver file, facets gene level file
    

    main:

    ch_versions = Channel.empty()

    ch_samplesheet.map {
            meta, maf, facets_gene, hla_file ->
                [meta, maf, hla_file]
                
        }
        .set { netMHCpan_input_ch }
    

    ch_samplesheet.map {
            meta, maf, facets_gene, hla_file ->
                [meta, maf, facets_gene]
                
        }
        .set { phylowgs_input_ch }

    // phylowgs workflow
    PHYLOWGS(phylowgs_input_ch)

    NETMHCPAN(netMHCpan_input_ch)

    phylowgs_output_ch = PHYLOWGS.out.summ.join(PHYLOWGS.out.muts, by:[0]).join(PHYLOWGS.out.mutass, by:[0])
    phylowgs_output_ch.view()


    NEOANTIGENINPUT(netMHCpan_input_ch,phylowgs_output_ch,NETMHCPAN.out.netmhcpanoutput)
    
    // NEOANTIGEN_EDITING()
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }



    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    neoin          = NEOANTIGENINPUT.out.json
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
