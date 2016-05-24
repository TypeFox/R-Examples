filename <- "rev.mod"
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, filename)
m <- matrix(nrow=3, ncol=3)
m[1,] <- c(1,2,3)
m[2,] <- c(1,5,10)
m[3,] <- c(10,4,2)
eq.freq <- c(1,2,3)
h <- hmm(m, eq.freq)
mod <- read.tm(filename)
mod2 <- mod
mod2$backgd <- rep(0.25, 4)
mod3 <- mod
mod3$backgd <- c(0.6, 0.1, 0.2, 0.1)
m <- simulate.msa(mod, 20)
m <- simulate.msa(list(mod, mod2, mod3), 20, hmm=h)
m <- matrix(1, nrow=3, ncol=3)
h <- hmm(m)
l <- simulate.msa(list(mod, mod2, mod3), 100, get.features=TRUE, hmm=h)
names(l)
l$msa
l$feats
coverage.feat(l$feats[l$feats$feature=="state1",])
unlink(filename)
