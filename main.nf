#!/usr/bin/env nextflow

nextflow.preview.dsl=2


file( "${params.outdir}" ).mkdirs()

println( "\n***\nParameters in use:")
params.each { println "${it}" }
println( "***\n")


/*
 * import modules
 */
include filter from './modules/filter' params(  //params)
    loom_input: {(params.containsKey('loom_input')) ? params.loom_input:'unfiltered.loom'},
    pyscenic_container: params.pyscenic_container ,
    outdir: params.outdir ,
    loom_filtered: 'filtered.loom' ,
    thr_min_genes: 200,
    thr_min_cells: 3,
    thr_n_genes: 5001,
    thr_pct_mito: 0.25,
)

// include './modules/pyscenic' params(params)
// include 'modules/BPanalysis'
// include 'modules/integrate'
    //(params.loom_output): 'pyscenic_integrated-output.loom'

/*
 * main workflow
 */
workflow {

    loomuf = Channel.fromPath( params.loom_input )
    loomf = filter( loomuf )

    //SCENIC(
    //    loomf[0],
    //    params.TFs
    //    )
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}


