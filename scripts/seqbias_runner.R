library(seqbias)
library(Rsamtools)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
args

ref_file <- args[1]
accession_file <- args[2] # i.e. MAY_A_1_accessions.txt
path <- args[3] # i.e. nobackup/data/MAY_A_1/
rm(args)

# for accession in accession_file
#   create reads, save, and plot files
#   run seqbias
#   

# Open reference file. Only need to do this once so outside loop
ref_fn <- system.file(ref_file)
ref_f <- FaFile(ref_file)
open.FaFile(ref_f)

accessions <- readLines(accession_file)

for (i in 1:length(accessions)) {
  accession = accessions[i]
  print(accession)


  # sample path: nobackup/data/MAY_A_1/{accession}/{accession}Aligned.out.sam.bam.sorted
  reads_file <- paste(path, "/", accession, "/", accession, "Aligned.out.sam.bam.sorted.bam", sep="")
  save_file <- paste(path, "/", accession, "/", accession, "-model.yml", sep="")
  plot_file <- paste(path, "/", accession, "/", accession, "-plot.png", sep="")
  corrected_plot_file <- paste(path, "/", accession, "/", accession, "-plotcorrected.png", sep="")
  biased_freqs_file <- paste(path, "/", accession, "/", accession, "-freqs.txt", sep="")
  corrected_freqs_file <- paste(path, "/", accession, "/", accession, "-corrected.txt", sep="")

  ## Assess bias of the sample
  
  # Generate random sample of intervals across ref genome
  print("generating random intervals...")
  ref_seqs <- scanFaIndex(ref_f)
  I <- random.intervals(ref_seqs, n=5000, ms=100000)
  print("...done")  

  # Extract sequences for these intervals from a FASTA file
  print("extracting sequences for these intervals from the reference")
  seqs <- scanFa(ref_f, I)
  neg_idx <- as.logical(I@strand == '-')
  seqs[neg_idx] <- reverseComplement(seqs[neg_idx])
  print("...done")  

  # Extract read counts across these intervals from a BAM file
  print("extracting read counts across these intervals")
  counts <- count.reads(reads_file, I, binary=T)
  print("...done")

  # Compute & plot kmer freqs
  print("computing kmer frequencies")
  freqs <- kmer.freq(seqs, counts)
  print("...done")
  print("writing kmer frequencies to file")
  write.table(freqs, biased_freqs_file)
  #lapply(freqs, write, biased_freqs_file, sep="\n")
  print("...done")
  P <- qplot(x = pos,
            y = freq,
            ylim = c(0.15, 0.4),
            color = seq,
            data = freqs,
            geom = "line")
  P <- P + facet_grid(seq ~ . )
  ggsave(filename=plot_file, width=12, height=8, units="cm")
  
  
  
  # Prediction
  print(ref_file)
  print(reads_file)
  sb <- seqbias.fit(ref_file, reads_file, L = 5, R = 15)
  print("finished fitting the model")
  seqbias.save(sb, save_file)
  print("finished saving")
  # sb <- seqbias.load(ref_file, save_file)
  bias <- seqbias.predict(sb, I)
  print("finished predicting")
  counts.adj <- mapply(FUN=`/`, counts, bias, SIMPLIFY = F)
  freqs.adj <- kmer.freq(seqs, counts.adj)
  write.table(freqs.adj, corrected_freqs_file)
  # lapply(freqs.adj, write, corrected_freqs_file, sep="\n")

  
  P <- qplot(x = pos,
            y = freq,
            ylim = c(0.15, 0.4),
            color = seq,
            data = freqs.adj,
            geom = "line")
  P <- P + facet_grid(seq ~ . )
  ggsave(filename=corrected_plot_file, width=12, height=8, units="cm")

}


