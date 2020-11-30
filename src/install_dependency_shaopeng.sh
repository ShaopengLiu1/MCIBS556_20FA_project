# to install necessary dependencies for the scripts
# prerequisites: conda or miniconda version >= 4.6
#
# Checklist:
#

#!/bin/bash
date
echo "Install all necessary dependencies for the project"



### local variables
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# to activate conda env
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
. ${conda_path}/etc/profile.d/conda.sh



### create conda env
conda create -y -p ${pipe_path}/MCIBS566_project_env_py37 python=3.7
conda activate ${pipe_path}/MCIBS566_project_env_py37
# install sourmash
conda install -y -c bioconda sourmash
# install kmer
conda install -y -c bioconda khmer
# install kmc
conda install -y -c bioconda kmc
# install sra toolkits
conda install -y -c bioconda sra-tools
# conda Bracken is NOT an official version, and is not properly installed by conda install
# install kracken2 in conda (bracken by git)
conda install -y -c bioconda kraken2

conda deactivate














echo "Finished!"
date
