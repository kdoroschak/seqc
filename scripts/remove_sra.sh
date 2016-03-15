#!/bin/bash
args=("$@")
path=${args[0]}
accessionsfile=${args[1]}

# set -e
while read p; do
		if [ -s $path/${p}_2.fastq ]; then 
			rm $path/${p}.sra
		fi
done < $accessionsfile

