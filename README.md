# Reproducible and scalable regulatory network reconstruction from single cell transcriptomics datasets using SCENIC

## Overview

## Requirements

The following tools are required to run the steps in this pipeline:
* [Nextflow](https://www.nextflow.io/)
* A container system, either of:
    * [Docker](https://docs.docker.com/)
    * [Singularity](https://www.sylabs.io/singularity/)

The following container images will be pulled by nextflow as needed:
* Docker: [aertslab/pyscenic:latest](https://hub.docker.com/r/aertslab/pyscenic).
* Singularity: [aertslab/pySCENIC:latest](https://www.singularity-hub.org/collections/2033).
* [See also here.](https://github.com/aertslab/pySCENIC#docker-and-singularity-images)


## Quick start

## Pipeline

### Starting point


### Filtering

Basic preprocessing filtering on both gene- and cell-level is done using tools included in the
[Scanpy](https://github.com/theislab/scanpy)
toolkit, which can efficiently deal with large datasets.
For this pipeline, we have re-packaged this software into a container with an associated Nextflow workflow that produces a filtered dataset with customizable thresholds.


### Analysis

#### Basic visualization with Scanpy

#### GRN Inference and construction of regulons

Gene regulatory network inference and regulon construction is performed using a python implementation of the SCENIC software package (pySCENIC).
A [Nextflow implementation of the SCENIC pipeline](https://github.com/aertslab/SCENICprotocol)
is used for this step.

Download/update the SCENIC repository:

    nextflow pull aertslab/SCENICprotocol
 
Run the SCENIC pipeline:

    nextflow run aertslab/SCENICprotocol \
        -profile singularity \
        --loom_input expr_mat.loom \
        --TFs allTFs_hg38.txt \
        --motifs motifs.tbl \
        --db *feather


### Visualization

#### Adding dimensionality reductions


#### SCope
SCope is a fast visualization tool for large-scale and high dimensional single-cell data.
Loom files produced by this pipeline can be viewed using the 
[SCope viewer](http://scope.aertslab.org).

To produce a loom file that is compatible with the SCope viewer:


## References and more information

### SCENIC
* [SCENIC Nextflow pipeline](https://github.com/aertslab/scenic-nf)
* [SCENIC (R) on GitHub](https://github.com/aertslab/SCENIC)
* [SCENIC website](http://scenic.aertslab.org/)
* [SCENIC publication](https://doi.org/10.1016/j.cell.2018.05.057)
* [pySCENIC on GitHub](https://github.com/aertslab/pySCENIC)
* [pySCENIC documentation](https://pyscenic.readthedocs.io/en/latest/)

### SCope
* [SCope webserver](http://scope.aertslab.org/)
* [SCope on GitHub](https://github.com/aertslab/SCope)
* [SCopeLoomR](https://github.com/aertslab/SCopeLoomR)
* [SCopeLoomPy](https://github.com/aertslab/SCopeLoomPy)

### Scanpy
* [Scanpy on GitHub](https://github.com/theislab/scanpy)
* [Scanpy documentation](https://scanpy.readthedocs.io/)
* [Scanpy publication](https://doi.org/10.1186/s13059-017-1382-0)




