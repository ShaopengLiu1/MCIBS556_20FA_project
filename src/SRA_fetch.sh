# range is the number of SRA assension numbers that represent individual samples.
# adjust the range number according to the length of the SRA list
# SRA_16s_list.txt, len = 260
# SRA_metagenomics_list. txt, len = 50

#!/bin/bash

for i in {1..50}
do
SRA=$(sed -n $( echo $i )p SRA_metagenomics_list.txt)
fastq-dump $SRA
SRA=0
done
