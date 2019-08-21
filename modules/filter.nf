
params.loom_input = '' // unfiltered loom input
params.loom_filtered = 'filtered.loom'
params.loom_output = 'pyscenic_integrated-output.loom'
params.thr_min_genes = 200
params.thr_min_cells = 3
params.thr_n_genes = 5000
params.thr_pct_mito = 0.25

/*
 * basic filtering
 */

process filter {
    cache 'deep'
    container params.pyscenic_container

    input:
    // file loomUnfiltered from file( params.loom_input )
    file loomuf

    output:
    file 'loom_filtered.loom'
    file 'anndata.h5ad'
    //file params.loom_filtered // into expr
    //file 'anndata.h5ad' // into SCfilter

    """
    filtering-basic.py \
        --loom_input $loomuf \
        --loom_filtered loom_filtered.loom \
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
workflow FILTER {
    get:
        loomuf
    main:
        filter( loomuf )
    emit:
        loomf.out
}

