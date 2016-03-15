#!/bin/bash
args=("$@")
path=${args[0]}
files=$path/*.sra
for f in $files; do
	echo Splitting SRA archive $p into fastq
	fastq-dump.2.5.4 --split-files $f -O $path
done