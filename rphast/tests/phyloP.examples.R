library("rphast")

## phyloP.sph (This is too slow for testing)
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
files <- c("ENr334-100k.fa", "rev.mod")
unzip(exampleArchive, files)
tm <- read.tm("rev.mod")
t1 <- phyloP.sph(tm, prior.only=TRUE, nsites=20)
plot(t1$nsub, t1$prior.distrib)
msa <- sub.msa(read.msa("ENr334-100k.fa", offset=41405894), end=10000)
t2 <- phyloP.sph(tm, msa)
t2
unlink(files)

rm(list = ls())
gc()
