#install.packages("entropy")
library("entropy")
library("ggplot2")

convertFreqsToMatrix <- function(freqs) {
  freq_a = freqs[freqs$seq == "a", ]$freq
  freq_c = freqs[freqs$seq == "c", ]$freq
  freq_g = freqs[freqs$seq == "g", ]$freq
  freq_t = freqs[freqs$seq == "t", ]$freq
  mat <- rbind(freq_a, freq_c, freq_g, freq_t)
  return(mat)
}

calculateKL <- function(freqs1, freqs2) {
  if ((dim(freqs1)[1] != dim(freqs2)[1]) || (dim(freqs1)[2] != dim(freqs2)[2])) {
    stop("dimensions of frequency matrices do not match")
  }
  ncols = dim(freqs1)[2]
  kl = sapply(1:ncols, function(col) KL.plugin(freqs1[,col], freqs2[,col], unit="log2"))
  return(kl) 
}

readFreqs <- function(filename) {
  return(read.table(filename))
}

frameKLResult <- function(kl) {
  pos <- seq(-floor(length(kl)/2), floor(length(kl)/2), 1)
  klresult <- data.frame(pos=pos, kl=kl)
  return(klresult)
}
# 
# s1 <- readFreqs("test.txt")
# s1 <- convertFreqsToMatrix(s1)
# 
# freqmat <- convertFreqsToMatrix(freqs)
# test <- calculateKL(freqmat, s1)
# test


args <- commandArgs(trailingOnly = TRUE)
print(args)

file1 <- args[1] #"SRR898239-freqs.txt" #
file2 <- args[2] #"SRR898240-freqs.txt" #

freqs1 <- convertFreqsToMatrix(readFreqs(file1))
freqs2 <- convertFreqsToMatrix(readFreqs(file2))

klresult <- frameKLResult(calculateKL(freqs1, freqs2))

P <- qplot(x = klresult$pos,
           y = klresult$kl,
           ylim = c(0, 0.5),
           #color = seq,
           data = klresult,
           geom = "line")
ggsave(filename="testklplot.png", width=12, height=8, units="cm")
