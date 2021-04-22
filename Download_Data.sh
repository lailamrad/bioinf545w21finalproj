#!/bin/bash
# would like to fetch all fastq files associated with SRP SRR13169280
# SRR sample range from *80 to *89 (10 samples total)
# return date and time 
echo "Beggan" $(date)
#this for loop sequentially downloads the raw ziped fastq files into the current dir 
for (( i = 0; i <= 9; i++))
do
	prefetch SRR1316928$i
	fastq-dump -I ----gzip --split-files SRR1316928$i
done
# return date and time 
echo "Finished" $(date)

