FROM python:3.6.8-slim
LABEL maintainer="christopher.flerin@kuleuven.vib.be"

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="build-essential apt-utils \
        python3-dev libhdf5-dev libfreetype6-dev libtool \
        m4 autoconf automake patch bison flex libpng-dev libopenblas-dev \
        tcl-dev tk-dev libxml2-dev zlib1g-dev libffi-dev cmake" && \
    apt-get update && \
    apt-get install -y debconf locales && dpkg-reconfigure locales && \
    apt-get install -y $BUILDPKGS && \
    ### run time:
    apt-get install -y zlib1g hdf5-tools gfortran libgcc1 libstdc++ musl \
        libopenblas-base tcl tk libxml2 libffi6 less procps

# install dependencies:
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /tmp/requirements.txt

RUN pip install --no-cache-dir pyscenic==0.9.9 && \
    pip install --no-cache-dir scanpy==1.4.3

RUN apt-get remove --purge -y $BUILDPKGS && \
    rm -rf /var/lib/apt/lists/*


