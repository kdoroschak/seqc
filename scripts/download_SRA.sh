#!/bin/bash
args=("$@")
filename=${args[0]}
path=${args[1]}
echo "Starting download"
while read p; do
	echo $p
	wget -P $path ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP025/SRP025982/$p/$p.sra
done < $filename
