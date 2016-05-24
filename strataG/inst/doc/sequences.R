## ----echo = FALSE, message = FALSE---------------------------------------
options(digits = 2)
library(strataG)

## ------------------------------------------------------------------------
data(dolph.seqs)

i <- sample(1:10, 1)
j <- sample(1:10, 1)
x <- c(rep("n", i), dolph.seqs[[1]], rep("n", j))
x
x.trimmed <- trimNs(as.DNAbin(x))
as.character(as.list(x.trimmed))

## ------------------------------------------------------------------------
bf <- baseFreqs(dolph.seqs)
bf$site.freqs[, 1:8]
bf$base.freqs

## ------------------------------------------------------------------------
fs <- fixedSites(dolph.seqs)
fs[1:20]

vs <- variableSites(dolph.seqs)
vs

## ------------------------------------------------------------------------
fs <- fixedSites(dolph.seqs, bases = c("c", "t"))
fs[1:20]

vs <- variableSites(dolph.seqs, bases = c("c", "t"))
vs

## ------------------------------------------------------------------------
iupacCode(c("c", "t", "t", "c", "c"))
iupacCode(c("c", "t", "a", "c", "c"))
iupacCode(c("g", "t", "a", "c", "c"))

## ------------------------------------------------------------------------
validIupacCodes(c("c", "t", "t", "c", "c"))
validIupacCodes(c("c", "t", "a", "c", "c"))
validIupacCodes(c("g", "t", "a", "c", "c"))

## ------------------------------------------------------------------------
createConsensus(dolph.seqs)

## ------------------------------------------------------------------------
nucleotideDiversity(dolph.seqs)

## ------------------------------------------------------------------------
# create gtypes
data(dolph.seqs)
data(dolph.strata)
dloop.haps <- cbind(dLoop = dolph.strata$id)
rownames(dloop.haps) <- dolph.strata$id
strata.schemes <- dolph.strata[, c("broad", "fine")]
rownames(strata.schemes) <- dolph.strata$id
dloop <- new("gtypes", gen.data = dloop.haps, ploidy = 1,
             schemes = strata.schemes, sequences = dolph.seqs,
             strata = "fine")
dloop <- labelHaplotypes(dloop, "Hap.")$gtypes

# calculate divergence
nucleotideDivergence(dloop)

## ------------------------------------------------------------------------

fixedDifferences(dloop)

## ------------------------------------------------------------------------
x <- as.DNAbin(dolph.seqs)
mostDistantSequences(x, num.seqs = 5)

## ------------------------------------------------------------------------
mostRepresentativeSequences(x, num.seqs = 5)

