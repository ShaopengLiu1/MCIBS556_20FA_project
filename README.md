# MCIBS556_20FA_project
---

#### Notes:
1. use abundance matrix for downstream analysis
2. data sharing:
   1. Shaopeng currently run data in his lab server, waiting for feedback from ICS helpdesk regarding sharing large files (~30GB) in the ICS server. 
   2. All scripts will be pushed to the github repo  
3. please put all finished scripts in the **src** folder and all intermediate files into a separate folder for each task  
   
#### Check list:

- [ ] 1. Data preparation:

  - [ ] 1.1 Pipeline: raw data -> count matrix (Bracken, kallastio) -> abundance matrix -> cleaned matrix
  - [ ] 1.2 Download NCBI bacteria genome/protein database 
  - [ ] 1.3 Pipeline: pull some files -> generate pairwise similarity matrix for a given k (for iterative tree building purpose)

- [ ] 2. Implement ML methods and performance evaluation for cleaned matrix
  - [ ] 2.1 Pick methods: random forest, ???
  - [ ] 2.2 Evaluation matrix
  - [ ] 2.3 pipeline: start from the cleaned matrix

- [ ] 3. Implement different DRs as control

  - [ ] 3.1 Taxonomic-based DR
  - [ ] 3.2 Data-driven DR
  - [ ] 3.3 What other papers use

- [ ] 4. Build similarity-based tree

  - [ ] 4.1 Tree logic: iteratively divide groups into subgroups (call pipe 1.3)?
  - [ ] 4.2 Combine 3 to build similairty-based DR

- [ ] 5. Results comparison
  
  - [ ] 5.1 Incorporate 3.1-3.3, 4.2 with pipe 2.3, compare the results

- [ ] 6. Extension (possibly the true meet for the future)

  - [ ] Compare similarity-based tree VS phylogenetic tree VS data-driven tree

