#!/usr/bin/env nextflow

println( "\n***\nParameters in use:")
params.each { println "${it}" }
println( "***\n")


// process importData {
// 
// }

file( "${params.outdir}" ).mkdirs()


/*
 * basic filtering
 */

process filter {

    input:
    file loomUnfiltered from file( params.loom_input )

    output:
    file params.loom_filtered into expr

    """
    filtering-basic.py \
        --loom_input ${loomUnfiltered} \
        --loom_filtered ${params.loom_filtered} \
        --thr_min_genes ${params.thr_min_genes} \
        --thr_min_cells ${params.thr_min_cells} \
        --thr_n_genes ${params.thr_n_genes} \
        --thr_pct_mito ${params.thr_pct_mito}
    """
}
expr.last().collectFile(storeDir:params.outdir)

/*
 * end of basic filtering
 */


/*
 * preprocess, visualize, project, cluster processing steps
 */

process preprocess {
    cache 'deep'
    input:
    file params.loom_filtered from expr
    output:
    file '01_preprocessed.h5ad' into SCpreprocess
    """
    preprocess_visualize_project_scanpy.py \
        preprocess \
        --loom_filtered ${params.loom_filtered} \
        --ad_preprocessed 01_preprocessed.h5ad \
        --threads ${params.threads}
    """
}

process pca {
    input:
    file '01_preprocessed.h5ad' from SCpreprocess
    output:
    file '02_pca.h5ad' into SCpca
    """
    preprocess_visualize_project_scanpy.py \
        pca \
        --ad_pca 02_pca.h5ad \
        --threads ${params.threads}
    """
}

process visualize {
    input:
    file '02_pca.h5ad' from SCpca
    output:
    file '03_visualize.h5ad' into SCvisualize
    """
    preprocess_visualize_project_scanpy.py \
        visualize \
        --ad_visualize 03_visualize.h5ad \
        --threads ${params.threads}
    """
}

process cluster {
    input:
    file '03_visualize.h5ad' from SCvisualize
    output:
    file '03_visualize.h5ad' into SCcluster
    """
    preprocess_visualize_project_scanpy.py \
        cluster \
        --ad_cluster 04_cluster.h5ad \
        --threads ${params.threads}
    """
}


/*
 * End of preprocess, visualize, project, cluster processing steps
 */


/*
 * SCENIC steps
 */

/*
 * end of SCENIC steps
 */



/*
 * results integration
 */

/*
 * end of results integration
 */

