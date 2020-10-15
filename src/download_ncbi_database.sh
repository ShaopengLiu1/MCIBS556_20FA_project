# download NCBI bacteria genome/protein database
# For simplicity: only download completed records, omitting all partial / incomplete data.
# This step would take few hours and ~60GB storage!
# ref: ftp://ftp.ncbi.nlm.nih.gov/genomes/



#!/bin/bash
date
echo "Downloading NCBI database, it may take several hours and ~60GB space!!!"

### local variables
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [ -d ${pipe_path}/NCBI_database ]; then
	echo "Resource folder NCBI_database already exists, please double check!"
else
	mkdir -p ${pipe_path}/NCBI_database/genomes
	mkdir ${pipe_path}/NCBI_database/proteins
	cd ${pipe_path}/NCBI_database
	#download the annotation file
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
	#





