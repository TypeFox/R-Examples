## ----echo = FALSE, message = FALSE---------------------------------------
library(strataG)

## ------------------------------------------------------------------------
#--- create a diploid (microsatellite) gtypes object
data(dolph.msats)
data(dolph.strata)
strata.schemes <- dolph.strata[, c("broad", "fine")]
rownames(strata.schemes) <- dolph.strata$id
msats <- new("gtypes", gen.data = dolph.msats[, -1], ploidy = 2,
  ind.names = dolph.msats[, 1], schemes = strata.schemes)

#--- create a haploid sequence (mtDNA) gtypes object
data(dolph.seqs)
dloop.haps <- cbind(haps = dolph.strata$id)
rownames(dloop.haps) <- dolph.strata$id
dloop <- new("gtypes", gen.data = dloop.haps, ploidy = 1, 
  sequences = dolph.seqs)

## ------------------------------------------------------------------------
msats.fine <- new("gtypes", gen.data = dolph.msats[, -1], ploidy = 2,
  ind.names = dolph.msats[, 1], schemes = strata.schemes, strata = "fine")

## ------------------------------------------------------------------------
# create a single data.frame with the msat data and stratification
msats.merge <- merge(dolph.strata, dolph.msats, all.y = TRUE, description = date())
str(msats.merge)

# create the gtypes object
msats.fine <- df2gtypes(msats.merge, ploidy = 2, id.col = 1, strata.col = 3, loc.col = 5)

## ------------------------------------------------------------------------
data(dolph.seqs)

seq.df <- dolph.strata[ c("id", "broad", "id")]
colnames(seq.df)[3] <- "D-loop"
dl.g <- df2gtypes(seq.df, ploidy = 1, sequences = dolph.seqs)
dl.g

## ------------------------------------------------------------------------
dl.haps <- labelHaplotypes(dl.g)
dl.haps$gtype

## ------------------------------------------------------------------------
data(dolph.haps)

haps.g <- sequence2gtypes(dolph.haps)
haps.g

## ------------------------------------------------------------------------
# extract and name the stratification scheme
strata <- dolph.strata$fine
names(strata) <- dolph.strata$ids

# create the gtypes object
dloop.fine <- sequence2gtypes(dolph.seqs, strata, seq.names = "dLoop",
  description = "dLoop: fine-scale stratification")
dloop.fine

## ------------------------------------------------------------------------
library(adegenet)
# from example(df2genind)
df <- data.frame(locusA=c("11","11","12","32"),
                 locusB=c(NA,"34","55","15"),
                 locusC=c("22","22","21","22"))
row.names(df) <- .genlab("genotype",4)
obj <- df2genind(df, ploidy=2, ncode=1)
obj

# convert to gtypes
gi.g <- genind2gtypes(obj)
gi.g

## ------------------------------------------------------------------------
sub.msats <- msats.fine[sample(nInd(msats.fine), 10), , ]
sub.msats

## ------------------------------------------------------------------------
sub.msats <- sub.msats[, c("D11t", "EV37", "TV7"), ]
sub.msats

## ------------------------------------------------------------------------
sub.msats <- msats.fine[, c("Ttr11", "D11t"), "Coastal"]
sub.msats

## ------------------------------------------------------------------------
msats.smry <- summary(msats.fine)
str(msats.smry)

## ------------------------------------------------------------------------
summarizeLoci(msats.fine)

## ------------------------------------------------------------------------
# randomly stratify individuals to two populations
strata(msats) <- sample(c("Pop1", "Pop2"), nInd(msats), rep = TRUE)
summary(msats)$strata.smry

## ------------------------------------------------------------------------
# choose "broad" stratification scheme
msats <- stratify(msats, "broad")
summary(msats)$strata.smry

## ------------------------------------------------------------------------
schemes(dloop) <- strata.schemes

## ------------------------------------------------------------------------
stratify(dloop, "fine")

## ------------------------------------------------------------------------
# unstratify a random 10 samples
x <- strata(msats)
x[sample(indNames(msats), 10)] <- NA
strata(msats) <- x
summary(msats)$strata.smry

## ------------------------------------------------------------------------
msats <- stratify(msats, "fine")

# original
summary(msats)$strata.smry

# permuted
ran.msats <- permuteStrata(msats)
summary(ran.msats)$strata.smry

## ------------------------------------------------------------------------
gen.mat <- as.matrix(msats)
head(gen.mat)

## ------------------------------------------------------------------------
gen.mat <- as.matrix(msats, one.col = TRUE)
head(gen.mat)

