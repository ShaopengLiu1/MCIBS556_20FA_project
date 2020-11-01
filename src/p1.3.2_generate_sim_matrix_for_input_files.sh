# Generate JI matrix by sourmash for query against ref
# input: query/ref containing absolute paths of genome files, range is the k_size to use
# output: a JI matrix
# Programmer: Shaopeng
# Last update: 10/31/2020



date
### read parameters
while getopts q:r:c:s:t:ah opts
do
	case "$opts" in
		q) query="$OPTARG";;		# query file containing all abs path
		r) ref="$OPTARG";;		# ref file containing all abs path
		c) range="$OPTARG";;		# range of size to check, format: start-end-gap
		s) scale="$OPTARG";;		# the scaling factor to use, default 2000
		t) threads="$OPTARG";;		# number of threads to use, default 24
		a) abundance="yes";;		# enable abundance tracking, default is no
		h) echo "
Benchmarking for sourmash
Usage: bash <script> -q <query> -r <ref> -c <range>

Parameters:
query/ref (-q/-r): files containing absolute path of all input files, in the output matrix, query is the row and ref is the column
range (-c): range of k values to run for the similarity matrix, e.g. 10-60-5
scale (-s): scaling factor in sourmash (proportion of MH to use), default 2000
threads (-t): number of threads to use, default 24
abundance (-a): enable abundance tracking (cosine similarity instead of JI), default is no
"
exit;;
[?]) echo "use -h for help"
exit;;
esac
done

# check input
if [ -z "$query" ] || [ -z "$ref" ] || [ -z "$range" ]; then
	echo "Missing input parameter!!! We need query / ref / range"
	exit 1
fi
[ -z "$scale" ] && scale=2000
[ -z "$threads" ] && threads=24
[ -z "$abundance" ] && abun_indicator="" || abun_indicator="--track-abundance" # for abundance option in building signature
# read in input
query=$(readlink -f $query)
ref=$(readlink -f $ref)
# range adjustment
temp_range=`echo $range | awk -F"-" '{ if(($1==1)) print $1+1"-"$2"-"$3; else print $1"-"$2"-"$3}'`
  r_start=`echo $temp_range | cut -d"-" -f 1`
  r_end=`echo $temp_range | cut -d"-" -f 2`
  r_gap=`echo $temp_range | cut -d"-" -f 3`
  r_adj_start=$((r_start+(r_end-r_start)%r_gap))
temp_range=${r_adj_start}-${r_end}-${r_gap}
kmer_sets=$(seq -s ","  $r_adj_start $r_gap $r_end) # for sourmash usage
echo "The kmer sets would be built on ${kmer_sets}"



### local variables
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# to activate conda env
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
. ${conda_path}/etc/profile.d/conda.sh
conda activate ${pipe_path}/MCIBS566_project_env_py37
# time cmd
ltime="/usr/bin/time -av -o temp_runLog"
# create output folder
time_tag=`date +"%m_%d_%H-%M"`
mkdir output_${time_tag}
cd output_${time_tag}



### build sourmash signatures of all input files (for our analysis: query = ref, so only need one)
# the threads option works for 10x bam files only
${ltime} sourmash compute -k ${kmer_sets} \
	--scaled ${scale} -p ${threads} ${abun_indicator} \
	-o ref.sig $(cat ${query} ${ref} | sort | uniq -c | awk '{print $2}' | paste -s)
mv temp_runLog record_sourmash_compute_ref.log
# genetrate matrix
for k in $(seq $r_adj_start $r_gap $r_end); do
	${ltime} sourmash compare -k ${k} -p  ${threads}  --csv JI_matrix_of_input_k${k}.csv ref.sig
	mv temp_runLog record_sourmash_compare_k${k}.log
done







