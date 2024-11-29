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
include { NETMHCSTABANDPAN } from '../subworkflows/msk/netmhcstabandpan/main'
include { NETMHCPAN } from '../modules/msk/netmhcpan/main'
include { NEOANTIGENUTILS_NEOANTIGENINPUT } from '../modules/msk/neoantigenutils/neoantigeninput'
include { NEOANTIGEN_EDITING } from '../subworkflows/msk/neoantigen_editing'
include { NEOANTIGENUTILS_CONVERTANNOTJSON } from '../modules/msk/neoantigenutils/convertannotjson'

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

    ch_cds_and_cdna = Channel.value([file(params.cds), file(params.cdna)])

    ch_samplesheet.map {
            meta, maf, facets_hisens_cncf, hla_file ->
                [meta, maf, hla_file]

        }
        .set { netMHCpan_input_ch }


    ch_samplesheet.map {
            meta, maf, facets_hisens_cncf, hla_file ->
                [meta, maf, facets_hisens_cncf]

        }
        .set { phylowgs_input_ch }

    ch_samplesheet.map {
            meta, maf, facets_hisens_cncf, hla_file ->
                [meta, []]

        }
        .set { ch_sv_empty }

    // phylowgs workflow
    PHYLOWGS(phylowgs_input_ch)

    ch_versions = ch_versions.mix(PHYLOWGS.out.versions)

    NETMHCSTABANDPAN(netMHCpan_input_ch,ch_cds_and_cdna,ch_sv_empty)

    ch_versions = ch_versions.mix(NETMHCSTABANDPAN.out.versions)

    netMHCpanMut = NETMHCSTABANDPAN.out.tsv
                        .filter{ it[0].typeMut == true && it[0].fromStab == false }
    netMHCpanWT = NETMHCSTABANDPAN.out.tsv
                        .filter{ it[0].typeMut == false && it[0].fromStab == false }
    stabNetMHCpanMut = NETMHCSTABANDPAN.out.tsv
                        .filter{ it[0].typeMut == true && it[0].fromStab == true }
    stabnetMHCpanWT = NETMHCSTABANDPAN.out.tsv
                        .filter{ it[0].typeMut == false && it[0].fromStab == true }

    merged = merge_for_input_generation(netMHCpan_input_ch, PHYLOWGS.out.summ, PHYLOWGS.out.muts, PHYLOWGS.out.mutass, netMHCpanMut, netMHCpanWT)

    merged_netMHC_input = merged
            .map{
                new Tuple(it[0], it[1], [], it[2])
            }
    merged_phylo_output = merged
        .map{
            new Tuple(it[0], it[3], it[4], it[5])
        }
    merged_netmhc_tsv = merged
        .map{
            new Tuple(it[0], it[6], it[7])
        }

    NEOANTIGENUTILS_NEOANTIGENINPUT(merged_netMHC_input,merged_phylo_output,merged_netmhc_tsv)

    ch_versions = ch_versions.mix(NEOANTIGENUTILS_NEOANTIGENINPUT.out.versions)

    NEOANTIGEN_EDITING(NEOANTIGENUTILS_NEOANTIGENINPUT.out.json, file(params.iedbfasta))

    ch_versions = ch_versions.mix(NEOANTIGEN_EDITING.out.versions)

    NEOANTIGENUTILS_CONVERTANNOTJSON(NEOANTIGEN_EDITING.out.annotated_output)

    ch_versions = ch_versions.mix(NEOANTIGENUTILS_CONVERTANNOTJSON.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }



    emit:
    versions         = ch_versions                 // channel: [ path(versions.yml) ]
    neo_out          = NEOANTIGEN_EDITING.out.annotated_output
    tsv_out          = NEOANTIGENUTILS_CONVERTANNOTJSON.out.neoantigenTSV
}

def merge_for_input_generation(netMHCpan_input_ch, summ_ch, muts_ch, mutass_ch, netmhcpan_mut_tsv_ch, netmhcpan_wt_tsv_ch ) {
    netMHCpan_input = netMHCpan_input_ch
        .map{
            new Tuple(it[0].id,it)
            }
    summ = summ_ch
        .map{
            new Tuple(it[0].id,it)
            }
    muts = muts_ch
        .map{
            new Tuple(it[0].id,it)
            }
    mutass = mutass_ch
        .map{
            new Tuple(it[0].id,it)
            }
    netmhcpan_mut_tsv = netmhcpan_mut_tsv_ch
        .map{
            new Tuple(it[0].id,it)
            }
    netmhcpan_wt_tsv = netmhcpan_wt_tsv_ch
        .map{
            new Tuple(it[0].id,it)
            }
    merged = netMHCpan_input
                .join(summ)
                .join(muts)
                .join(mutass)
                .join(netmhcpan_mut_tsv)
                .join(netmhcpan_wt_tsv)
                .map{
                    new Tuple(it[1][0], it[1][1], it[1][2], it[2][1], it[3][1], it[4][1], it[5][1], it[6][1])
                }
    return merged
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

