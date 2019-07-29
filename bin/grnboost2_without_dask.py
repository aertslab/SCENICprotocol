#!/usr/bin/env python3

import sys
import time
import loompy as lp
import pandas as pd
from multiprocessing import Pool, cpu_count
import argparse

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2, _prepare_input
from arboreto.core import SGBM_KWARGS, RF_KWARGS, EARLY_STOP_WINDOW_LENGTH
from arboreto.core import to_tf_matrix, target_gene_indices, infer_partial_network


################################################################################
################################################################################

parser_grn = argparse.ArgumentParser(description='Run GRNBoost2 using a multiprocessing pool')

parser_grn.add_argument('expression_mtx_fname',
        type=argparse.FileType('r'),
        help='The name of the file that contains the expression matrix for the single cell experiment.'
        ' Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).')
parser_grn.add_argument('tfs_fname',
        type=argparse.FileType('r'),
        help='The name of the file that contains the list of transcription factors (TXT; one TF per line).')
parser_grn.add_argument('-o', '--output',
        type=argparse.FileType('w'), default=sys.stdout,
        help='Output file/stream, i.e. a table of TF-target genes (CSV).')
parser_grn.add_argument('-t', '--transpose', action='store_const', const = 'yes',
        help='Transpose the expression matrix (rows=genes x columns=cells).')
parser_grn.add_argument('--num_workers',
        type=int, default=cpu_count(),
        help='The number of workers to use. (default: {}).'.format(cpu_count()))
parser_grn.add_argument('--seed', type=int, required=False, default=None,
        help='Seed value for regressor random state initialization (optional)')

parser_grn.add_argument('--cell_id_attribute',
        type=str, default='CellID',
        help='The name of the column attribute that specifies the identifiers of the cells in the loom file.')
parser_grn.add_argument('--gene_attribute',
        type=str, default='Gene',
        help='The name of the row attribute that specifies the gene symbols in the loom file.')

args = parser_grn.parse_args()


################################################################################
################################################################################
################################################################################


def runInferPartialNet( target_gene_index ):

    target_gene_name = gene_names[target_gene_index]
    target_gene_expression = expression_matrix[:, target_gene_index]

    n = infer_partial_network(
        regressor_type='GBM',
        regressor_kwargs=SGBM_KWARGS,
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        target_gene_name=target_gene_name,
        target_gene_expression=target_gene_expression,
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=args.seed)
    return( n )

if __name__ == '__main__':

    lf = lp.connect( args.expression_mtx_fname.name, mode='r', validate=False )
    # genes in columns:
    ex_matrix = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID ).T
    lf.close()
    gene_names = ex_matrix.columns
    print('Loaded expression matrix of {} cells and {} genes...'.format( ex_matrix.shape[0], ex_matrix.shape[1] ) , file=sys.stderr)

    tf_names = load_tf_names( args.tfs_fname.name )
    print('Loaded {} TFs...'.format( len(tf_names) ) , file=sys.stderr)

    expression_matrix, gene_names, tf_names = _prepare_input(ex_matrix, gene_names, tf_names)
    tf_matrix, tf_matrix_gene_names = to_tf_matrix(expression_matrix, gene_names, tf_names)

    print('starting GRNBoost2 using {} processes...'.format( args.num_workers ), file=sys.stderr)
    start_time = time.time()

    with Pool( args.num_workers ) as p:
        adjs = p.map(runInferPartialNet, 
                target_gene_indices(gene_names, target_genes='all')
                )
    adj = pd.concat( adjs ).sort_values(by='importance', ascending=False)

    end_time = time.time()
    print('Done in {} seconds.'.format(end_time - start_time), file=sys.stderr)
    adj.to_csv( args.output, index=False, sep="\t")

