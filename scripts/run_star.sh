#!/bin/bash
args=("$@")
path=${args[0]}
accessionsfile=${args[1]}

# set -e
while read p; do
	echo $p #accession is $p
	mkdir $path/${p}
	#echo Splitting SRA archive into fastq
	#fastq-dump.2.5.4 --split-files -O $path $path/${p}.sra || exit 1
	#rm $path/$p.sra
	STAR --genomeDir ~/star-genome --readFilesIn $path/${p}_1.fastq $path/${p}_2.fastq --runThreadN 30 --outFileNamePrefix $path/${p}/$p --outSAMtype BAM SortedByCoordinate || exit 1
	./software/tophat-2.1.0.Linux_x86_64/samtools_0.1.18 index $path/$p/${p}Aligned.sortedByCoord.out.bam
done < $accessionsfile
Rscript freq_runner.R GRCh38.p3_genomic.fa $accessionsfile $path
