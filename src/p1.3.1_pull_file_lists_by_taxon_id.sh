# generate abs file paths for input taxon_id 
# input: a txt file in which each line represent a taxon_id
# output: a txt file containing abs path of all corresponding genome files
# programmer: Shaopeng
# last update: 11/29/2020



date
input_file=$1
echo "The input file is ${input_file}"
[ -z $input_file ] && echo "Missing input taxon list!!!" && exit 1
### remove trailing blank lines if any
sed -i '/^$/d' ${input_file}

### specific operation for "otutable.txt" input in this analysis
if [ $input_file == "otutable.txt" ]; then
	cut -f 1 otutable.txt | sed '1d' | cut -d"_" -f 4 > species_list.txt
	input_file=species_list.txt
fi


echo "Note: the associated-species file uses species taxon id, which contains sub-string information"
echo "Only complete genomes would be used, all contigs would be ignored"
### considerations:
# 1. the full NCBI database contains 757440 records with 21273 full genomes, the 1st size is too big to download
# 2. contigs have various size and multiple submissions: it's very hard to pick an appropriate one; and the size would directly affect JI (small overlap) and leads to unexpected artifacts



### local variables
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
record_file=${pipe_path}/NCBI_database/NCBI_complete_genome_quick_list.txt
total=$(wc -l $input_file | awk '{print $1}')
count=0
full=$(cat ${input_file} | wc -l)


### grep record
# when encounter dups: 1) use representative genome if any; 2) use newer submission if none representative. That's why use a sort -k2,2rV there
grep -w -f ${input_file} ${record_file} | sort -k2,2rV  > full_grep_record_with_dups.txt
echo -n > genome_file_path.txt
echo -n > protein_file_path.txt

for spe_id in $(cat ${input_file}); do
	found=$(grep ${spe_id} full_grep_record_with_dups.txt | wc -l)
	if (( $found > 0 )); then
		let "count+=1"
		accession=$(grep ${spe_id} full_grep_record_with_dups.txt | head -1 | cut -f 1)
		f_name=$(ls ${pipe_path}/NCBI_database/genomes | grep $accession)
		echo "${pipe_path}/NCBI_database/genomes/${f_name}" >> genome_file_path.txt
		f2_name=$(ls ${pipe_path}/NCBI_database/proteins | grep $accession)
		echo "${pipe_path}/NCBI_database/proteins/${f2_name}" >> protein_file_path.txt
	fi
done
# rm empty record (those with no matched files but the folder path only) in protain db
sed -i '/NCBI_database\/proteins\/$/d' protein_file_path.txt
sed -i '/NCBI_database\/genomes\/$/d' genome_file_path.txt #this shouldn't happen, but just in case



### check missing rate
missing=$(python -c "print(1-${count}*1.0/${total})")
echo "In total, ${count} out of ${total} genomes have been found in the complete genome list."
echo "The missing rate is ${missing}"



### end
echo "pipeline finished"
date









