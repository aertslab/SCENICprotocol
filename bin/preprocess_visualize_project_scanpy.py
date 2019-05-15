#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import argparse

################################################################################
################################################################################

def preprocess( args ):
    sc.logging.print_versions()
    # read filtered data from a loom file
    adata = sc.read_loom(args.loom_filtered)

    # Total-count normalize (library-size correct) to 10,000 reads/cell
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

    # log transform the data.
    sc.pp.log1p(adata)

    # identify highly variable genes.
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata)

    # keep only highly variable genes:
    adata = adata[:, adata.var['highly_variable']]

    # mito and genes/counts cuts
    mito_genes = adata.var_names.str.startswith('MT-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1

    # regress out total counts per cell and the percentage of mitochondrial genes expressed
    sc.pp.regress_out(adata, ['n_counts', 'percent_mito'], n_jobs=args.threads)

    # scale each gene to unit variance, clip values exceeding SD 10.
    sc.pp.scale(adata, max_value=10)

    adata.write( args.anndata )


def pca( args ):
    adata = sc.read_h5ad( args.anndata )
    # principal component analysis
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True)
    adata.write( args.anndata )

def pcaNdimSelect( args ):
    print( "pcaNdim" )

def visualize( args ):
    adata = sc.read_h5ad( args.anndata )
    # neighborhood graph of cells (determine optimal number of PCs here)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    # compute UMAP
    sc.tl.umap(adata)
    adata.write( args.anndata )


def cluster( args ):
    adata = sc.read_h5ad( args.anndata )
    # cluster the neighbourhood graph
    sc.tl.louvain(adata,resolution=0.4)

    sc.pl.umap(adata, color=['louvain'] )
    # find marker genes
    sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

    # sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
    # sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
    adata.write( args.anndata )


def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.as_matrix()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


def integrateOutput( args ):
    print( "output" )

################################################################################
################################################################################

FUNCTIONS = {
    'preprocess': preprocess ,
    'pca': pca ,
    'pcaNdimSelect': pcaNdimSelect ,
    'visualize': visualize ,
    'cluster': cluster ,
    'integrateOutput': integrateOutput ,
    }

parser = argparse.ArgumentParser(description='Preprocess, visualize, project using Scanpy.')
parser.add_argument( 'command', choices=FUNCTIONS.keys() )
parser.add_argument('--loom_filtered', help='Loom file with basic filtering', required=False, default='filtered.loom' )
parser.add_argument('--threads', help='Maximum number of threads to use', required=False, default=6, type=int )

parser.add_argument('--anndata', help='Intermediate filename storing Scanpy preprocessing output', required=False, default='anndata.h5ad' )
args = parser.parse_args()
func = FUNCTIONS[args.command]

################################################################################
################################################################################

if __name__ == "__main__":
    func( args )
    # preprocess( args )
    # pca( adataf="01_preprocessed.h5ad", args )
    # visualize( adataf="02_pca.h5ad", args )
    # cluster( adataf="03_visualize.h5ad", args )
    # integrateOutput( adataf="04_cluster.h5ad", args )

