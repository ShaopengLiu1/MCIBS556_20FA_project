# MCIBS556_20FA_project
To repeat all the results presented for MCBIS556 project, please follow the steps below:  
## Install dependencies
1. Install conda or [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (version >= 4.6)  
2. Clone this repo  
```
git clone https://github.com/ShaopengLiu1/MCIBS556_20FA_project.git
```
3. Go to `MCIBS556_20FA_project/src` folder, then run the scripts listed below:  
- install conda environment by (!!!this file will be merged with Morris' and be renamed later!!!)
```
bash install_dependency_shaopeng.sh
``` 
- download the NCBI bacteira genome/protein database 
```
#Note: this script would run ~10 hours and need ~40GB storage!!!
#For convenience, "nohup" is recomended, please run the cmd and just leave it to run (it's safe to quit your session with nohup)
nohup bash download_ncbi_database.sh &
```



