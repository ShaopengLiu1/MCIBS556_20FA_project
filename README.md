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

## Repeat results in each steps
1. Build abundance matrix from input NGS data in fastq format
2. To be continued......


## Project notes
1. use abundance matrix for downstream analysis
2. data sharing:
   1. Shaopeng currently run data in his lab server
   2. How to share files in ICS: email both user id to the i-ask center and the admin will put the 2 accounts into one group
   3. All scripts will be pushed to the github repo  
3. please put all finished scripts in the **src** folder and all intermediate files into a separate folder for each task  
   
## Progress checklist
- [ ] 1. Resources  

  - [x] 1.1 Download NCBI bacteria genome/protein database 
  - [ ] 1.2 Find bacteria lists for specific targets (biological traits, human microbiome project, whole phylogenetic tree)

- [ ] 2. Prepare input data (by week of 10.26)  

  - [ ] 2.1 Pipeline: produce abundance matrix from input metagenome/16S data
  - [ ] 2.2 Pipeline: generate a similarity matrix for a given list of bacteria

- [ ] 3. Implement different DRs as control 

  - [ ] 3.1 Pipeline: Taxonomic-based DR
  - [ ] 3.2 Pipeline: Data-driven DR based on covariance
  - [ ] 3.3 What other papers use (need to check)

- [ ] 4. Build similarity-based tree  

  - [ ] 4.1 Tree logic: iteratively divide groups into subgroups
  - [ ] 4.2 Combine 3 to build similairty-based DR
  
- [ ] 5. Implement ML methods and performance evaluation for cleaned matrix 

  - [ ] 5.1 Pick methods: random forest, neuron network, GLM, etc.
  - [ ] 5.2 Evaluation matrix
  - [ ] 5.3 pipeline: run abundance matrix with some DR method and generate the evaluation matrix

- [ ] 6. Results comparison  
  
  - [ ] 6.1 Incorporate 3.1-3.3, 4.2 with pipe 5.3, compare the results
  - [ ] 6.2 Writing!

- [ ] 7. Extension (possibly the true meet for the future)  

  - [ ] Compare similarity-based tree VS phylogenetic tree VS data-driven tree

