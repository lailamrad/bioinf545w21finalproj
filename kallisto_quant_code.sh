#!/bin/bash

# running ettertivly not sure if this is the best thoguh
for ID in `tail -n2 sampIDlist22.txt`
do
	kallisto quant -i ./mus_musculus/transcriptome.idx -o ./kallisto_out/SRR${ID}_kallisto_out --single -l 60 -s 6 ./Raw_data/SRR${ID}_*.fastq.gz
done
echo "--------script is done-----------"

