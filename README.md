# MCIBS556_20FA_project
To repeat all the results presented for MCBIS556 project, please follow the steps below:  
## Install dependencies
1. Install conda or [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (version >= 4.6)  
2. Clone this repo  
```
git clone https://github.com/ShaopengLiu1/MCIBS556_20FA_project.git
```
3. Go to `MCIBS556_20FA_project/src` folder, then run the scripts listed below:
   3.1. install conda environment by (!!!this file will be merged with Morris' and be renamed later!!!)
   ```
   bash install_dependency_shaopeng.sh
   ``` 
   3.2. download the NCBI bacteira genome/protein database 
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
- [ ] 1. Data preparation (week of 10.12)

  - [ ] 1.1 Pipeline: raw data -> count matrix (Bracken, kallastio) -> abundance matrix -> cleaned matrix
  - [x] 1.2 Download NCBI bacteria genome/protein database 
  - [ ] 1.3 Pipeline: pull some files -> generate pairwise similarity matrix for a given k (for iterative tree building purpose)

- [ ] 2. Implement different DRs as control (week of 10.19 and 10.26)

  - [ ] 2.1 Taxonomic-based DR
  - [ ] 2.2 Data-driven DR
  - [ ] 2.3 What other papers use

- [ ] 3. Implement ML methods and performance evaluation for cleaned matrix (week of 11.2)
  - [ ] 3.1 Pick methods: random forest, neuron network, GLM, etc.
  - [ ] 3.2 Evaluation matrix
  - [ ] 3.3 pipeline: start from the cleaned matrix

- [ ] 4. Build similarity-based tree (week of 11.9 and 11.16)

  - [ ] 4.1 Tree logic: iteratively divide groups into subgroups (call pipe 1.3)?
  - [ ] 4.2 Combine 3 to build similairty-based DR

- [ ] 5. Results comparison (week of 11.23 and 11.30)
  
  - [ ] 5.1 Incorporate 3.1-3.3, 4.2 with pipe 2.3, compare the results
  - [ ] 5.2 Writing!

- [ ] 6. Extension (possibly the true meet for the future)

  - [ ] Compare similarity-based tree VS phylogenetic tree VS data-driven tree

