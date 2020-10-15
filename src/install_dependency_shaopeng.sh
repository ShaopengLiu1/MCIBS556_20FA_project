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

### create conda env
conda create -y -p ${pipe_path}/MCIBS566_project_env_py37 python=3.7



echo "Finished!"
date
