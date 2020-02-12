#!/usr/bin/env nextflow

println( "\n***\nParameters in use:")
params.each { println "${it}" }
println( "***\n")


file( "${params.outdir}" ).mkdirs()


/*
 * basic filtering
 */

process filter {
    cache 'deep'
    container params.pyscenic_container

    input:
    file loomUnfiltered from file( params.loom_input )

    output:
    file params.loom_filtered into expr
    file 'anndata.h5ad' into SCfilter

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
    container params.pyscenic_container
    input:
    file 'anndata.h5ad' from SCfilter
    output:
    file 'anndata.h5ad' into SCpreprocess
    """
    preprocess_visualize_project_scanpy.py \
        preprocess \
        --loom_filtered ${params.loom_filtered} \
        --anndata anndata.h5ad \
        --threads ${params.threads}
    """
}

process pca {
    cache 'deep'
    container params.pyscenic_container
    input:
    file 'anndata.h5ad' from SCpreprocess
    output:
    file 'anndata.h5ad' into SCpca
    """
    preprocess_visualize_project_scanpy.py \
        pca \
        --anndata anndata.h5ad \
        --threads ${params.threads}
    """
}

process visualize {
    cache 'deep'
    container params.pyscenic_container
    input:
    file 'anndata.h5ad' from SCpca
    output:
    file 'anndata.h5ad' into SCvisualize
    """
    preprocess_visualize_project_scanpy.py \
        visualize \
        --anndata anndata.h5ad \
        --threads ${params.threads}
    """
}

process cluster {
    cache 'deep'
    container params.pyscenic_container
    input:
    file 'anndata.h5ad' from SCvisualize
    output:
    file 'anndata.h5ad' into SCcluster
    """
    preprocess_visualize_project_scanpy.py \
        cluster \
        --anndata anndata.h5ad \
        --threads ${params.threads}
    """
}

/*
 * End of preprocess, visualize, project, cluster processing steps
 */


/*
 * SCENIC steps
 */

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

process GRNinference {
    cache 'deep'
    container params.pyscenic_container

    input:
    file TFs from tfs
    file params.loom_filtered from expr
    // file exprMat from expr

    output:
    file 'adj.tsv' into GRN

    """
    arboreto_with_multiprocessing.py \
        ${params.loom_filtered} \
        ${TFs} \
        --method ${params.grn} \
        --num_workers ${params.threads} \
        -o adj.tsv \
        ${(params.containsKey('sparse')) ? '--sparse' : ''} \
        ${(params.containsKey('seed')) ? "--seed  ${params.seed}" : ""} \
        --cell_id_attribute ${params.cell_id_attribute} \
        --gene_attribute ${params.gene_attribute} \
    """
}

process cisTarget {
    cache 'deep'
    container params.pyscenic_container

    input:
    file exprMat from expr
    file 'adj.tsv' from GRN
    file feather from featherDB
    file motif from motifs

    output:
    file 'reg.csv' into regulons

    """
    pyscenic ctx \
        adj.tsv \
        ${feather} \
        --annotations_fname ${motif} \
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

    input:
    file exprMat from expr
    file 'reg.csv' from regulons

    output:
    file params.pyscenic_output into scenicAUC1, scenicAUC2

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

/*
 * end of SCENIC steps
 */



/*
 * results integration
 */

process integrateOutput {
    cache 'deep'
    container params.pyscenic_container

    input:
    file params.pyscenic_output from scenicAUC2
    file 'anndata.h5ad' from SCcluster
    file 'scenic_umap.txt' from aucDRumap
    file 'scenic_tsne.txt' from aucDRtsne

    output:
    file params.loom_output into finalLoom

    """
    integrateOutput.py \
        --anndata anndata.h5ad \
        --loom_pyscenic ${params.pyscenic_output} \
        --loom_output ${params.loom_output} \
    """
}
finalLoom.last().collectFile(storeDir:params.outdir)

/*
 * end of results integration
 */

