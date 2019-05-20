#!/usr/bin/env python3

import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import argparse

import zlib
import json
import base64
import datetime

################################################################################
################################################################################

parser = argparse.ArgumentParser(description='Basic filtering using Scanpy')
parser.add_argument('--anndata', help='Intermediate filename storing Scanpy preprocessing output', required=True, default='anndata.h5ad' )
parser.add_argument('--loom_pyscenic', help='Loom file from pySCENIC', required=True, default='pyscenic.loom' )
parser.add_argument('--loom_output', help='Final loom file with pySCENIC and Scanpy results integrated', required=True, default='pyscenic.loom' )
args = parser.parse_args()

################################################################################
################################################################################

def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.as_matrix()]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr

def integrateOutput( args ):
    adata = sc.read_h5ad( args.anndata )
    lf = lp.connect( args.loom_pyscenic, mode='r', validate=False )

    tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])

    Embeddings_X = pd.DataFrame()
    Embeddings_Y = pd.DataFrame()

    Embeddings_X["1"] = pd.DataFrame(adata.obsm['X_umap'])[0]
    Embeddings_Y["1"] = pd.DataFrame(adata.obsm['X_umap'])[1]

    Embeddings_X["2"] = pd.DataFrame(adata.obsm['X_pca'])[0]
    Embeddings_Y["2"] = pd.DataFrame(adata.obsm['X_pca'])[1]

    pc_to_use = 2

    metaJson = {}

    metaJson['embeddings'] = [
        {
            "id": -1,
            "name": f"Scanpy t-SNE {pc_to_use}PC"
        },
        {
            "id": 1,
            "name": f"Scanpy UMAP {pc_to_use}PC"
        },
        {
            "id": 2,
            "name": "Scanpy PC1/PC2"
        }
    ]

    metaJson["clusterings"] = [{
                "id": 0,
                "group": "Scanpy",
                "name": "Scanpy louvain default resolution",
                "clusters": [],
            }]

    metaJson["metrics"] = [
            {
                "name": "nUMI"
            }, {
                "name": "nGene"
            }, {
                "name": "Percent_mito"
            }
    ]

    metaJson["annotations"] = [
        #{
        #    "name": "Genotype",
        #    "values": list(set(adata.obs['Genotype'].values))
        #},
        #{
        #    "name": "Timepoint",
        #    "values": list(set(adata.obs['Timepoint'].values))
        #},
        #{
        #    "name": "Sample",
        #    "values": list(set(adata.obs['Sample'].values))
        #}
    ]


    for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
        clustDict = {}
        clustDict['id'] = i
        clustDict['description'] = f'Unannotated Cluster {i + 1}'
        metaJson['clusterings'][0]['clusters'].append(clustDict)

    clusterings = pd.DataFrame()

    clusterings["0"] = adata.obs['louvain'].values.astype(np.int64)

    col_attrs = {
        "CellID": np.array(adata.obs.index),
        "nUMI": np.array(adata.obs['n_counts'].values),
        "nGene": np.array(adata.obs['n_genes'].values),
        #"Genotype": np.array(adata.obs['Genotype'].values),
        #"Timepoint": np.array(adata.obs['Timepoint'].values),
        #"Sample": np.array(adata.obs['Sample'].values),
        "Percent_mito": np.array(adata.obs['percent_mito'].values),
        "Embedding": dfToNamedMatrix(tsneDF),
        "Embeddings_X": dfToNamedMatrix(Embeddings_X),
        "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
        "Clusterings": dfToNamedMatrix(clusterings),
        "ClusterID": np.array(adata.obs['louvain'].values)
    }

    row_attrs = {
        "Gene": np.array(adata.var.index)
    }

    attrs = {
        "title": "sampleTitle",
        "MetaData": json.dumps(metaJson),
        "Genome": 'hg38',
        "SCopeTreeL1": "",
        "SCopeTreeL2": "",
        "SCopeTreeL3": ""
    }

    attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

    #f'/staging/leuven/stg_00003/cbd-bioinf/{project}/{sampleName}.loom'
    #adata.raw.X
    #datetime.date.today().strftime('%Y%m%d')

    lp.create(
        filename = args.loom_output ,
        #layers=pd.DataFrame(adata.X.todense()).T.values, 
        layers=pd.DataFrame(adata.X).T.values, 
        row_attrs=row_attrs, 
        col_attrs=col_attrs, 
        file_attrs=attrs
    )



if __name__ == "__main__":
    integrateOutput( args )

