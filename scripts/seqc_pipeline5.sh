#!/bin/bash
args=("$@")
path=${args[0]}
accessionsfile=${args[1]}

function maxthreads {
	while [ `jobs | wc -l` -ge 6 ]
	do
		sleep 30
	done
}

process_accession ()
{
	accession=$1

	echo Downloading accession ${accession}.
	wget -nv -P $path ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP025/SRP025982/${accession}/${accession}.sra

	echo Splitting fastq for ${accession}.
	fastq-dump.2.5.4 --split-files ${path}/${accession}.sra -O $path
	
        if [ -s $path/${accession}_2.fastq ]; then
                echo Splitting fastq  was successful for accession ${accession}, so removing SRA files.
                rm $path/${accession}.sra
        fi

	echo Running STAR for accession ${accession}.
	mkdir ${path}/${accession}
	STAR --genomeDir ~/star-genome --readFilesIn $path/${accession}_1.fastq $path/${accession}_2.fastq --runThreadN 30 --outFileNamePrefix $path/${accession}/${accession} --outSAMtype BAM SortedByCoordinate
        
	echo Indexing bam for accession ${accession}.
	./software/tophat-2.1.0.Linux_x86_64/samtools_0.1.18 index $path/${accession}/${accession}Aligned.sortedByCoord.out.bam	
	if [ -s $path/${accession}/${accession}Aligned.sortedByCoord.out.bam ]; then
		echo STAR was successful for accession ${accession}, so removing fastq files.
		rm $path/${accession}_1.fastq
		rm $path/${accession}_2.fastq
	fi
	
	echo Done processing accession ${accession}.
}

while read accession; do
        maxthreads
	echo --- Beginning to process accession ${accession}.
        process_accession $accession &
done < $accessionsfile
wait

echo ---------- Done processing all files, calculating frequencies. ----------
Rscript freq_runner.R GRCh38.p3_genomic.fa $accessionsfile $path

echo Done!!
