#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import argparse

################################################################################
################################################################################

parser = argparse.ArgumentParser(description='Basic filtering using Scanpy')
parser.add_argument('--loom_input', help='Unfiltered loom file', required=True, default='unfiltered.loom' )
parser.add_argument('--loom_filtered', help='Filtered loom file', required=True, default='filtered.loom' )
parser.add_argument('--anndata', help='Intermediate filename storing Scanpy preprocessing output', required=False, default='anndata.h5ad' )
parser.add_argument('--thr_min_genes', help='Threshold on minimum genes expressed per cell', required=False, default=200, type=int )
parser.add_argument('--thr_min_cells', help='Threshold on minimum cells in which a gene is detected', required=False, default=3, type=int )
parser.add_argument('--thr_n_genes', help='Threshold on maximum total counts per cell', required=False, default=5000, type=int )
parser.add_argument('--thr_pct_mito', help='Threshold on maximum fraction of counts in mito genes vs all genes', required=False, default=0.25, type=float )
args = parser.parse_args()

################################################################################
################################################################################


def initialFiltering( args ):

    sc.logging.print_versions()

    ##################################################
    # read unfiltered data from a loom file
    adata = sc.read_loom(args.loom_input)

    ##################################################
    # basic filtering / stats

    nCountsPerGene = np.sum(adata.X, axis=0)
    nCellsPerGene = np.sum(adata.X>0, axis=0)

    # Show info
    print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - " ,nCountsPerGene.max())
    print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - " ,nCellsPerGene.max())

    nCells=adata.X.shape[0]

    # pySCENIC thresholds
    minCountsPerGene=3*.01*nCells # 3 counts in 1% of cells
    print("minCountsPerGene: ", minCountsPerGene)

    minSamples=.01*nCells # 1% of cells
    print("minSamples: ", minSamples)

    ####################
    # initial cuts
    sc.pp.filter_cells(adata, min_genes=args.thr_min_genes)
    sc.pp.filter_genes(adata, min_cells=args.thr_min_cells)

    ####################
    # mito and genes/counts cuts
    mito_genes = adata.var_names.str.startswith('MT-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    if( sum(mito_genes)==0 ):
        adata.obs['percent_mito'] = 0.0
    else:
        adata.obs['percent_mito'] = np.ravel(np.sum(np.asmatrix(adata[:, mito_genes].X.todense()), axis=1)) / np.ravel(np.sum(adata.X, axis=1))
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))

    adata = adata[adata.obs['n_genes'] < args.thr_n_genes, :]
    adata = adata[adata.obs['percent_mito'] < args.thr_pct_mito, :]

    ##################################################
    # output to loom file:
    row_attrs = {
        "Gene": np.array(adata.var_names) ,
    }
    col_attrs = {
        "CellID": np.array(adata.obs_names) ,
        "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
        "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
    }

    lp.create(args.loom_filtered, adata.X.transpose(), row_attrs, col_attrs)
    adata.write( args.anndata )

if __name__ == "__main__":
    initialFiltering( args )

