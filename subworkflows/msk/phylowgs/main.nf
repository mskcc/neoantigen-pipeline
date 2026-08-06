include { PHYLOWGS_CREATEINPUT } from '../../../modules/msk/phylowgs/createinput/main'
include { PHYLOWGS_PARSECNVS } from '../../../modules/msk/phylowgs/parsecnvs/main'
include { PHYLOWGS_MULTIEVOLVE } from '../../../modules/msk/phylowgs/multievolve/main'
include { PHYLOWGS_WRITERESULTS } from '../../../modules/msk/phylowgs/writeresults/main'
include { PHYLOWGS_MULTIEVOLVEGO } from '../../../modules/local/phylowgs/multievolvego/main'
include { PHYLOWGS_WRITERESULTSGO } from '../../../modules/local/phylowgs/writeresultsgo/main'

workflow PHYLOWGS {

    take:
    ch_input_maf_and_genelevel

    main:

    ch_versions = channel.empty()

    ch_genelevel = ch_input_maf_and_genelevel
                        .map{
                            [it[0],it[2]]
                        }

    ch_maf = ch_input_maf_and_genelevel
                        .map{
                            [it[0],it[1]]
                        }

    PHYLOWGS_PARSECNVS(ch_genelevel)

    ch_versions = ch_versions.mix(PHYLOWGS_PARSECNVS.out.versions)

    ch_maf_and_cnv = join_maf_with_cnv(ch_maf,PHYLOWGS_PARSECNVS.out.cnv)

    PHYLOWGS_CREATEINPUT(ch_maf_and_cnv)

    ch_versions = ch_versions.mix(PHYLOWGS_CREATEINPUT.out.versions)

    // params.phylowgs_engine: 'python' (default) | 'go' | 'both'
    //   python -> run only the Python/C++ multievolve (original behavior)
    //   go     -> run only the Go reimplementation; its output becomes summ/muts/mutass
    //   both   -> run both; summ/muts/mutass stay Python, go_* carries the Go side for diffing
    def run_python = params.phylowgs_engine in ['python', 'both']
    def run_go     = params.phylowgs_engine in ['go', 'both']

    ch_summ      = channel.empty()
    ch_muts      = channel.empty()
    ch_mutass    = channel.empty()
    ch_go_summ   = channel.empty()
    ch_go_muts   = channel.empty()
    ch_go_mutass = channel.empty()

    if (run_python) {
        PHYLOWGS_MULTIEVOLVE(PHYLOWGS_CREATEINPUT.out.phylowgsinput)

        ch_versions = ch_versions.mix(PHYLOWGS_MULTIEVOLVE.out.versions)

        PHYLOWGS_WRITERESULTS(PHYLOWGS_MULTIEVOLVE.out.trees)

        ch_versions = ch_versions.mix(PHYLOWGS_WRITERESULTS.out.versions)

        ch_summ   = PHYLOWGS_WRITERESULTS.out.summ
        ch_muts   = PHYLOWGS_WRITERESULTS.out.muts
        ch_mutass = PHYLOWGS_WRITERESULTS.out.mutass
    }

    if (run_go) {
        PHYLOWGS_MULTIEVOLVEGO(PHYLOWGS_CREATEINPUT.out.phylowgsinput)

        ch_versions = ch_versions.mix(PHYLOWGS_MULTIEVOLVEGO.out.versions)

        PHYLOWGS_WRITERESULTSGO(PHYLOWGS_MULTIEVOLVEGO.out.trees)

        ch_versions = ch_versions.mix(PHYLOWGS_WRITERESULTSGO.out.versions)

        ch_go_summ   = PHYLOWGS_WRITERESULTSGO.out.summ
        ch_go_muts   = PHYLOWGS_WRITERESULTSGO.out.muts
        ch_go_mutass = PHYLOWGS_WRITERESULTSGO.out.mutass

        if (!run_python) {
            // go-only: Go results are the primary output feeding the rest of the pipeline
            ch_summ   = ch_go_summ
            ch_muts   = ch_go_muts
            ch_mutass = ch_go_mutass
        }
    }

    emit:

    summ        = ch_summ                               // channel: [ val(meta), [ summ ] ]
    muts        = ch_muts                                // channel: [ val(meta), [ muts ] ]
    mutass      = ch_mutass                              // channel: [ val(meta), [ mutass ] ]
    go_summ     = ch_go_summ                            // channel: [ val(meta), [ summ ] ]   (empty unless phylowgs_engine includes go)
    go_muts     = ch_go_muts                            // channel: [ val(meta), [ muts ] ]   (empty unless phylowgs_engine includes go)
    go_mutass   = ch_go_mutass                          // channel: [ val(meta), [ mutass ] ] (empty unless phylowgs_engine includes go)
    versions    = ch_versions                           // channel: [ versions.yml ]
}

def join_maf_with_cnv(maf,cnv) {
        def maf_channel = maf
            .map{
                [it[0].id,it]
                }
        def cnv_channel = cnv
            .map{
                [it[0].id,it]
                }
        def mergedWithKey = maf_channel
            .join(cnv_channel)
        def merged = mergedWithKey
            .map{
                [it[1][0],it[1][1],it[2][1]]
            }
        return merged

}
