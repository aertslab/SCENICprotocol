
# Installation of pySCENIC and dependencies

## Container system
Our preference is to use a container system to run pySCENIC, as it's more simple to install and update.
However, the drawbacks are that it's then difficult to include additional packages without re-building the entire image (note that `pip install --user ...` can sometimes be used to circumvent this).
See also this [additional information](https://github.com/aertslab/pySCENIC#docker-and-singularity-images).

### Docker

A pre-built container is available from DockerHub:
    [aertslab/pyscenic:latest](https://hub.docker.com/r/aertslab/pyscenic):
```bash
docker pull aertslab/pyscenic:0.9.18
```

This image can also be built from scratch on your own system using the pySCENIC Dockerfile and requirements file.
The pySCENIC version needs to be specified during the build:
```bash
wget https://raw.githubusercontent.com/aertslab/pySCENIC/master/Dockerfile
wget https://raw.githubusercontent.com/aertslab/pySCENIC/master/requirements_docker.txt
docker build -t aertslab/pyscenic:0.9.18 . --build-arg version=0.9.18
```
The build time is estimated to be approximately 10 minutes on a modern desktop.

Basic test of the image:
```bash
docker run -it --rm aertslab/pyscenic:0.9.18 pyscenic -h
```

### Singularity
Singularity images are currently hosted on [Singularity Hub](https://singularity-hub.org).
These can be obtained by:
```bash
singularity pull --name aertslab-pyscenic-0.9.18.sif shub://aertslab/pySCENIC:0.9.18
```
    
However, Singularity Hub currently requires an account for most actions, including container pulls.
Alternatives are:

#### Build the Singularity container from the Docker image
* From the public DockerHub:
```bash
singularity build aertslab-pyscenic-0.9.18.sif docker://aertslab/pyscenic:0.9.18
```

* From the local Docker daemon (this will only work if you have already built or pulled the docker image to your local machine):
```bash
singularity build aertslab-pyscenic-0.9.18.sif docker-daemon://aertslab/pyscenic:0.9.18
```
The build time is estimated to be approximately 10 minutes on a modern desktop.

#### Build the container from the [recipe file](https://github.com/aertslab/pySCENIC/blob/master/Singularity):
```bash
wget https://raw.githubusercontent.com/aertslab/pySCENIC/master/Singularity.0.9.18
wget https://raw.githubusercontent.com/aertslab/pySCENIC/master/requirements_docker.txt
singularity build aertslab-pyscenic-0.9.18.sif Singularity.0.9.18
```

Basic test of the image:
```bash
singularity run aertslab-pyscenic-0.9.18.sif pyscenic -h
```

### Volume mounts
By default, later versions of Singularity automatically mount the user home folder inside of the container and any files within are accessible to the container software.
However, when using Docker, or when needing to access files outside of the home folder on Singularity, it is necessary to manually specify the volumes to be mounted.
The easiest method is to specify the same common path on the host and container.
For example, if trying to access two resource files at:
```bash
/data/resources/cisTarget/motifs/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
/data/resources/cisTarget/databases/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
```
Mounting the `/data` volume is sufficient.
The container would need to be started as:
```bash
singularity run -B /data:/data aertslab-pyscenic-0.9.18.sif pyscenic [...]
```
or:
```bash
docker run -it --rm -v /data:/data aertslab/pyscenic:0.9.18 pyscenic [...]
```

## Conda

Conda can also be used to easily install the necessary packages:

```bash
conda create -n scenic_protocol python=3.6
conda activate scenic_protocol
```

Install some basic dependencies:

```bash
conda install numpy pandas matplotlib seaborn
conda install -c anaconda cytoolz
```

Install [Scanpy](https://scanpy.readthedocs.io/en/latest/installation.html):

```bash
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph louvain
conda install -c conda-forge multicore-tsne
pip install scanpy
```

Install pySCENIC

```bash
pip install pyscenic
```

Install environment as kernel for Jupyter

```bash
pip install --user ipykernel
python -m ipykernel install --user --name=scenic_protocol
```

Basic test of the environment:
```bash
pyscenic -h
```

## Interactive use with a Jupyter notebook

Both the container and the conda environments can be used as a remote kernel in a Jupyter notebook.

For example, to add a Singulariy kernel to Jupyter Lab installation, run a command similar to the following from the Jupyter Lab server:
```bash
python3 -m remote_ikernel manage --add \
    --kernel_cmd="singularity run -B /data /path/to/aertslab-pyscenic-0.9.18.sif ipython kernel -f {connection_file}" \
    --name="pyscenic-0918-singularity" \
    --interface=ssh \
    --host=hostname \
    --workdir="~/" \
    --language=python3
```
Note that the `-B` option needs to be used here to specify which files from the host system will be available inside the contaienr.

A similar command allows the use of a Conda environment in place of Singularity:
```bash
python3 -m remote_ikernel manage --add \
    --kernel_cmd="source /path/to/anaconda3/etc/profile.d/conda.sh && conda activate scenic_protocol && /path/to/envs/scenic_protocol/bin/ipython3 kernel -f {connection_file}" \
    --name="pyscenic-0918-conda" \
    --interface=ssh \
    --host=hostname \
    --workdir="~/" \
    --language=python3
```

After this, it should be possible to select the pySCENIC kernel in a Jupyter notebook.

