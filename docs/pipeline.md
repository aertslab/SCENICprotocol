
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

