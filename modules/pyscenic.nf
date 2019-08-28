
include './filter' params(params)

/*
 * SCENIC parameters
 */

// pySCENIC resources:
// transcription factors
params.TFs = "allTFs_hg38.txt"

// cisTarget resources (motif-based)
// Motif-to-TF annotation:
params.motifANN = "resources/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
// cisTarget ranking database
params.motifDB = "resources/hg19*mc9nr.feather"

// cisTarget resources (based on epigenomic tracks). Optional addition
// Track-to-TF annotation:
params.trackANN = ""
// cisTarget ranking database
params.trackDB = ""

params.grn = "grnboost2"
// params.pyscenic_output = "pyscenic_output.loom"
params.loom_filtered = 'filtered.loom'

// parameters for reading loom files:
params.cell_id_attribute = "CellID"
params.gene_attribute = "Gene"

type = '' // labels regulon, AUC matrix as motif or track

/*
 * SCENIC steps
 */

/*
// channel for SCENIC databases resources:
featherDB = Channel
    .fromPath( params.db )
    .collect() // use all files together in the ctx command

n = Channel.fromPath(params.db).count().get()
if( n==1 ) {
    println( "***\nWARNING: only using a single feather database:\n  ${featherDB.get()[0]}.\nTo include all database files using pattern matching, make sure the value for the '--db' parameter is enclosed in quotes!\n***\n" )
} else {
    println( "***\nUsing $n feather databases:")
    featherDB.get().each {
        println "  ${it}"
    }
    println( "***\n")
}

// expr = file(params.expr)
tfs = file(params.TFs)
motifs = file(params.motifs)
*/

process GRNinference_Dask {
    cache 'deep'
    container params.pyscenic_container

    input:
    file params.loom_filtered
    file params.TFs

    output:
    file 'adj.tsv'

    """
    pyscenic grn \
        --num_workers ${params.threads} \
        -o adj.tsv \
        --method ${params.grn} \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute} \
        ${params.loom_filtered} \
        ${params.TFs}
    """
}

process GRNinference_woDask {
    cache 'deep'
    container params.pyscenic_container

    input:
    file params.loom_filtered
    file params.TFs

    output:
    file 'adj.tsv'

    """
    grnboost2_without_dask.py \
        --output adj.tsv \
        --num_workers ${params.threads} \
        ${params.loom_filtered} \
        ${params.TFs}
    """
}

process cisTarget {
    cache 'deep'
    container params.pyscenic_container

    input:
    file exprMat
    file 'adj.tsv'
    file feather
    file ann
    val type

    output:
    file 'reg_${type}.csv'

    """
    pyscenic ctx \
        adj.tsv \
        ${feather} \
        --annotations_fname ${ann} \
        --expression_mtx_fname ${exprMat} \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute} \
        --mode "dask_multiprocessing" \
        --output reg.csv \
        --num_workers ${params.threads} \
    """
}

process AUCell {
    cache 'deep'
    container params.pyscenic_container
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file exprMat
    file 'reg.csv'

    output:
    file params.pyscenic_output

    """
    pyscenic aucell \
        $exprMat \
        reg.csv \
        -o ${params.pyscenic_output} \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute} \
        --num_workers ${params.threads}
    """
}


/*
process visualizeAUC {
    cache 'deep'
    container params.pyscenic_container

    input:
    file params.pyscenic_output from scenicAUC1

    output:
    file 'scenic_umap.txt' into aucDRumap
    file 'scenic_tsne.txt' into aucDRtsne

    """
    preprocess_visualize_project_scanpy.py \
        visualizeAUC \
        --loom_pyscenic ${params.pyscenic_output} \
    """
}
*/

/*
 * end of SCENIC steps
 */

workflow SCENIC { //( loomf, tfs ) {
    get:
        loomf
        tfs
    main:
        /* GRN */
        tfs = file(params.TFs)
        //grn = GRNinference_Dask( loomf, tfs )
        grn = GRNinference_woDask( loomf, tfs )

        /* cisTarget */
        // channel for SCENIC databases resources:
        motifDB = Channel
            .fromPath( params.db )
            .collect() // use all files together in the ctx command
        motifs = file(params.motifs)

        grn.view()
        ctx_mtf = cisTarget( loomf, grn, params.motifDB, params.motifANN, 'mtf' )
        //ctx_trk = cisTarget( loomf, grn, trackDB, trackANN, 'trk' )

        /* AUCell */
        auc_mtf = AUCell( loomf, ctx_mtf )
        auc_mtf.view()
        //auc_trk = AUCell( loomf, ctx_trk )

        //visualizeAUC

}

