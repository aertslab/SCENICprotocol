
//thr_min_genes = 200
//thr_min_cells = 3
//thr_n_genes = 5000
//thr_pct_mito = 0.25

/*
 * basic filtering
 */

process filter {
    cache 'deep'
    container params.pyscenic_container
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file loomuf

    output:
    file params.loom_filtered
    file 'anndata.h5ad'

    """
    filtering-basic.py \
        --loom_input $loomuf \
        --loom_filtered ${params.loom_filtered} \
        --thr_min_genes ${params.thr_min_genes} \
        --thr_min_cells ${params.thr_min_cells} \
        --thr_n_genes ${params.thr_n_genes} \
        --thr_pct_mito ${params.thr_pct_mito}
    """
}

/*
 * end of basic filtering
 */


// workflow
// workflow FILTER {
//     get:
//         loomuf
//     main:
//         filter( loomuf )
//     emit:
//         loomf.out
// }

