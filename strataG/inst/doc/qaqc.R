## ----echo = FALSE, message = FALSE---------------------------------------
options(digits = 2)
library(strataG)

## ------------------------------------------------------------------------
data(msats.g)
msats <- msats.g
smry <- summarizeLoci(msats)
head(smry)

## ------------------------------------------------------------------------
# Find samples that share alleles at 2/3rds of the loci
dupGenotypes(msats, num.shared = 0.66)

## ------------------------------------------------------------------------
data(dolph.seqs)
seq.smry <- summarizeSeqs(as.DNAbin(dolph.seqs))
head(seq.smry)

## ------------------------------------------------------------------------
bf <- baseFreqs(as.DNAbin(dolph.seqs))

# nucleotide frequencies by site
bf$site.freq[, 1:15]

# overall nucleotide frequencies
bf$base.freqs

## ------------------------------------------------------------------------
lowFreqSubs(as.DNAbin(dolph.seqs), min.freq = 2)

## ------------------------------------------------------------------------
data(dolph.haps)
haplotypeLikelihoods(as.DNAbin(dolph.haps))

