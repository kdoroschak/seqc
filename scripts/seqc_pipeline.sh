#!/bin/bash
# Run as ./seqc_pipeline.sh <output directory> <accessions file>
#  ex: ./seqc_pipeline.sh ~/nobackup/AGR_A_1 A_accessions/AGR_A_1_accessions.txt

# Other notes:
# I installed samtools at: ./software/tophat-2.1.0.Linux_x86_64/samtools_0.1.1
#   It seems like I must have had to get the binaries through the tophat 
#   distribution (can't remember). Let me know if you have issues, but looking
#   up this particular distribution of tophat should be straightforward.
#
# Save to the ~/nobackup directory so the systems people don't get too mad

args=("$@")
path=${args[0]}
accessionsfile=${args[1]}

# Run at most 5 jobs at once.
# My bottleneck here was how many aligners to run at once.
# There are smarter ways to do this, but it doesn't impact the time too much.  
function maxjobs {
	while [ `jobs | wc -l` -ge 5 ]
	do
		sleep 30
	done
}

# process_accession does the following for a single NCBI SRA accession:
#  - Download the .sra file 
#  - Split the .sra into paired-end fastq
#  - If the split was successful, get rid of the .sra files
#  - Align 
#  - Index the alignment
#  - If the alignment was successful, get rid of the paired-end fastq  
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

	# Replace your aligner here
	echo Running STAR for accession ${accession}.
	mkdir ${path}/${accession}
	STAR --genomeDir ~/star-genome --readFilesIn $path/${accession}_1.fastq $path/${accession}_2.fastq --runThreadN 30 --outFileNamePrefix $path/${accession}/${accession} --outSAMtype BAM SortedByCoordinate
        
	echo Indexing bam for accession ${accession}.
	./software/tophat-2.1.0.Linux_x86_64/samtools_0.1.18 index $path/${accession}/${accession}Aligned.sortedByCoord.out.bam	

	# !!! Important !!! If you might be running another alignment on these
	# files, keep the fastq's around. 
	if [ -s $path/${accession}/${accession}Aligned.sortedByCoord.out.bam ]; then
		echo STAR was successful for accession ${accession}, so removing fastq files.
		rm $path/${accession}_1.fastq
		rm $path/${accession}_2.fastq
	fi
	
	echo Done processing accession ${accession}.
}

# Read in the file of accessions, and for each, do: 
while read accession; do
	# Check if we can run more jobs
        maxjobs
	# Process the next accession from the file.
	echo --- Beginning to process accession ${accession}.
        process_accession $accession &
done < $accessionsfile
wait

# Once all files are done, do any necessary summary calculations.
# "wait" command waits for all jobs spawned by the while loop to finish.
echo ---------- Done processing all files, calculating frequencies. ----------
Rscript freq_runner.R GRCh38.p3_genomic.fa $accessionsfile $path

echo Done!!
