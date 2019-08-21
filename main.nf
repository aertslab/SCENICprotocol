#!/usr/bin/env nextflow

nextflow.preview.dsl=2


loomuf = Channel.fromPath( params.loom_input )
file( "${params.outdir}" ).mkdirs()

println( "\n***\nParameters in use:")
params.each { println "${it}" }
println( "***\n")


/*
 * import modules
 */
include 'modules/filter' params(params)
// include 'modules/pyscenic'
// include 'modules/BPanalysis'
// include 'modules/integrate'

workflow {
    filter( loomuf )
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}


