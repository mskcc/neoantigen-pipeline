/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_neoqual_nf_pipeline'
include { PHYLOWGS_CREATEINPUT } from '../modules/msk/phylowgs/createinput/main'
include { PHYLOWGS_MULTIEVOLVE } from '../modules/msk/phylowgs/multievolve/main'
include { PHYLOWGS_PARSECNVS } from '../modules/msk/phylowgs/parsecnvs/main'
include { PHYLOWGS_WRITERESULTS } from '../modules/msk/phylowgs/writeresults/main'
include { PHYLOWGS } from '../subworkflows/msk/phylowgs'
include { GENERATE_MUTATED_PEPTIDES } from '../subworkflows/msk/generate_mutated_peptides/main'
include { NETMHCSTABANDPAN } from '../subworkflows/msk/netmhcstabandpan/main'
include { NEOANTIGENUTILS_NEOANTIGENINPUT } from '../modules/msk/neoantigenutils/neoantigeninput'
include { NEOANTIGEN_EDITING } from '../subworkflows/msk/neoantigen_editing'
include { NEOANTIGENUTILS_CONVERTANNOTJSON } from '../modules/msk/neoantigenutils/convertannotjson'
include { PHYLOWGS_STUB } from '../modules/local/phylowgs/stub/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NEOQUAL {

    take:
    ch_samplesheet // channel: samplesheet read in from --input It should have maf, polysolver file, facets gene level file
    outdir


    main:

    ch_versions = Channel.empty()

    ch_gtf_and_cdna = Channel.value([file(params.gtf), file(params.cdna)])

    // Generate mutated peptides (MUTALYZER_RETRIEVER + GENERATEMUTFASTA + GENERATEHLASTRING + NEOSV)
    ch_maf_hla_sv = ch_samplesheet.map { meta, maf, facets_hisens_cncf, hla_file ->
        [meta, maf, hla_file, null]  // null sv = no structural variants
    }

    GENERATE_MUTATED_PEPTIDES(
        ch_maf_hla_sv,
        Channel.value(file(params.reference_fasta)),
        Channel.value(file(params.reference_gff3)),
        Channel.value(file(params.gtf)),
        Channel.value(file(params.cdna))
    )
    ch_versions = ch_versions.mix(GENERATE_MUTATED_PEPTIDES.out.versions)

    ch_fasta_and_hla = GENERATE_MUTATED_PEPTIDES.out.mut_fasta
        .join(GENERATE_MUTATED_PEPTIDES.out.wt_fasta)
        .join(GENERATE_MUTATED_PEPTIDES.out.hla_string)

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

    ch_sv = GENERATE_MUTATED_PEPTIDES.out.sv_mut_fasta
        .join(GENERATE_MUTATED_PEPTIDES.out.sv_wt_fasta)

    // phylowgs workflow (optional)
    if ( params.run_phylowgs ) {
        PHYLOWGS(phylowgs_input_ch)
        ch_versions = ch_versions.mix(PHYLOWGS.out.versions)
        phylowgs_summ   = PHYLOWGS.out.summ
        phylowgs_muts   = PHYLOWGS.out.muts
        phylowgs_mutass = PHYLOWGS.out.mutass
    } else {
        PHYLOWGS_STUB(phylowgs_input_ch)
        ch_versions = ch_versions.mix(PHYLOWGS_STUB.out.versions)
        phylowgs_summ   = PHYLOWGS_STUB.out.summ
        phylowgs_muts   = PHYLOWGS_STUB.out.muts
        phylowgs_mutass = PHYLOWGS_STUB.out.mutass
    }

    NETMHCSTABANDPAN(ch_fasta_and_hla, ch_sv)

    ch_versions = ch_versions.mix(NETMHCSTABANDPAN.out.versions)

    netMHCpanMut = NETMHCSTABANDPAN.out.tsv
                        .filter{ it[0].typeMut == true && it[0].fromStab == false }
    netMHCpanWT = NETMHCSTABANDPAN.out.tsv
                        .filter{ it[0].typeMut == false && it[0].fromStab == false }
    stabNetMHCpanMut = NETMHCSTABANDPAN.out.tsv
                        .filter{ it[0].typeMut == true && it[0].fromStab == true }
    stabnetMHCpanWT = NETMHCSTABANDPAN.out.tsv
                        .filter{ it[0].typeMut == false && it[0].fromStab == true }

    merged = merge_for_input_generation(netMHCpan_input_ch, phylowgs_summ, phylowgs_muts, phylowgs_mutass, netMHCpanMut, netMHCpanWT)

    merged_netMHC_input = merged
            .map{
                [it[0], it[1], [], it[2]]
            }
    merged_phylo_output = merged
        .map{
            [it[0], it[3], it[4], it[5]]
        }
    merged_netmhc_tsv = merged
        .map{
            [it[0], it[6], it[7]]
        }

    NEOANTIGENUTILS_NEOANTIGENINPUT(merged_netMHC_input,merged_phylo_output,merged_netmhc_tsv,ch_gtf_and_cdna)

    ch_versions = ch_versions.mix(NEOANTIGENUTILS_NEOANTIGENINPUT.out.versions)

    NEOANTIGEN_EDITING(NEOANTIGENUTILS_NEOANTIGENINPUT.out.json, file(params.iedbfasta))

    ch_versions = ch_versions.mix(NEOANTIGEN_EDITING.out.versions)

    NEOANTIGENUTILS_CONVERTANNOTJSON(NEOANTIGEN_EDITING.out.annotated_output)

    ch_versions = ch_versions.mix(NEOANTIGENUTILS_CONVERTANNOTJSON.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'neoqual_software_versions.yml',
            sort: true,
            newLine: true
        )


    emit:
    versions         = ch_versions                 // channel: [ path(versions.yml) ]
    neo_out          = NEOANTIGEN_EDITING.out.annotated_output
    tsv_out          = NEOANTIGENUTILS_CONVERTANNOTJSON.out.neoantigenTSV
}

def merge_for_input_generation(netMHCpan_input_ch, summ_ch, muts_ch, mutass_ch, netmhcpan_mut_tsv_ch, netmhcpan_wt_tsv_ch ) {
    def netMHCpan_input = netMHCpan_input_ch
        .map{
            [it[0].id,it]
            }
    def summ = summ_ch
        .map{
            [it[0].id,it]
            }
    def muts = muts_ch
        .map{
            [it[0].id,it]
            }
    def mutass = mutass_ch
        .map{
            [it[0].id,it]
            }
    def netmhcpan_mut_tsv = netmhcpan_mut_tsv_ch
        .map{
            [it[0].id,it]
            }
    def netmhcpan_wt_tsv = netmhcpan_wt_tsv_ch
        .map{
            [it[0].id,it]
            }
    def merged = netMHCpan_input
                .join(summ)
                .join(muts)
                .join(mutass)
                .join(netmhcpan_mut_tsv)
                .join(netmhcpan_wt_tsv)
                .map{
                    [it[1][0], it[1][1], it[1][2], it[2][1], it[3][1], it[4][1], it[5][1], it[6][1]]
                }
    return merged
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
