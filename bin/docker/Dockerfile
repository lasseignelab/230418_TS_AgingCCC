FROM bioconductor/bioconductor_docker:RELEASE_3_16

###----------------------------------------------------------------------------##
##  R PACKAGES
##----------------------------------------------------------------------------##

#### Leave the top line the same for each new project. Uncomment lines to add ####
RUN R -e 'install.packages(c("tidyverse", "remotes",  "BiocManager", "devtools", "knitr", "lintr", "styler", "here", "markdown", "rmarkdown"))' \
  && R -e 'install.packages(c("ggplot2", "circlize", "reshape", "intergraph", "pheatmap", "clustree", "fields", "enrichR", "UpSetR", "scCustomize", "hdf5r", "rliger", "NMF", "gprofiler2", "GGally"))' \
  && R -e 'devtools::install_github("saeyslab/nichenetr")' \
  && R -e 'devtools::install_github("saezlab/liana")' \
  && R -e 'devtools::install_github("eddelbuettel/harmony")' \
  && R -e 'devtools::install_github("satijalab/sctransform", ref = "develop")' \
  && R -e 'remotes::install_github(repo = "satijalab/seurat", ref = "develop")' \
  && R -e 'BiocManager::install(c("STRINGdb", "Nebulosa", "DropletUtils", "dittoSeq", "glmGamPoi", "SummarizedExperiment", "ComplexHeatmap", "limma"))'

##----------------------------------------------------------------------------##
## SYSTEM PACKAGES
##----------------------------------------------------------------------------##

RUN apt-get update --yes && \
  apt-get upgrade --yes

RUN apt-get install --yes \
  python3-pip \
  python3-dev \
  python3-pip \
  python3-tk


RUN apt-get install --yes \
  libcurl4-openssl-dev \
  libssl-dev \
  libcairo2-dev libxt-dev \
  libxml2-dev \
  libudunits2-dev \
  libhdf5-dev \
  libv8-dev \
  libgdal-dev \
  xorg libx11-dev libglu1-mesa-dev \
  samtools libboost-all-dev \
  libgsl-dev

RUN apt-get install -y \
        wget \
        unzip \
        git


##----------------------------------------------------------------------------##
## PYTHON PACKAGES
##----------------------------------------------------------------------------##
  

RUN pip3 install leidenalg

  
##----------------------------------------------------------------------------##
## 3RD PARTY PROGRAMS
##----------------------------------------------------------------------------##
  
#### Uncomment these lines if you need a conda env in your docker ####
  
# Install miniconda
# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
#     /bin/bash ~/miniconda.sh -b -p /opt/conda
# ENV PATH=$CONDA_DIR/bin:$PATH

        
##----------------------------------------------------------------------------##
## CONDA ENV
##----------------------------------------------------------------------------##

#### Uncomment these lines if you need a conda env in your docker ####

#RUN conda config --set channel_priority false && \
#    conda env create -f opt/BLAZE-1.1.0/conda_env/environment.yml

#SHELL ["conda", "run", "-n", "blaze", "/bin/bash", "-c"]
