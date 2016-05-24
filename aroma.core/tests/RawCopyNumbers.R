library("aroma.core")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number of loci
J <- 500

mu <- double(J)
mu[100:150] <- mu[100:150] + 1
mu[320:400] <- mu[320:400] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y)))


cn <- RawCopyNumbers(y, x)
print(cn)

cn2 <- extractSubset(cn, subset=xSeq(cn, by=5))
print(cn2)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot along genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(cn, ylim=c(-3,3))
title(main="Complete and subsetted loci")
points(cn2, col="red", pch=176, cex=2)

legend("topright", pch=c(19,176), col=c("#999999", "red"), sprintf(c("raw [n=%d]", "every 5th [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cn2))), bty="n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Binned smoothing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(cn, col="#999999", ylim=c(-3,3))
title(main="Binned smoothing")

cnSa <- binnedSmoothing(cn, by=3)
lines(cnSa, col="blue")
points(cnSa, col="blue")

cnSb <- binnedSmoothing(cn, by=9)
lines(cnSb, col="red")
points(cnSb, col="red")

legend("topright", pch=19, col=c("#999999", "blue", "red"), sprintf(c("raw [n=%d]", "Bin(w=3) [n=%d]", "Bin(w=9) [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cnSa), nbrOfLoci(cnSb))), bty="n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Binned smoothing with shuffled loci
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
idxs <- sample(1:nbrOfLoci(cn))
cnR <- cn[idxs,]
plot(cnR, col="#999999", ylim=c(-3,3))
title(main="Binned smoothing (shuffled loci)")

cnSa <- binnedSmoothing(cnR, by=3)
lines(cnSa, col="blue")
points(cnSa, col="blue")

cnSb <- binnedSmoothing(cnR, by=9)
lines(cnSb, col="red")
points(cnSb, col="red")

legend("topright", pch=19, col=c("#999999", "blue", "red"), sprintf(c("raw [n=%d]", "Bin(w=3) [n=%d]", "Bin(w=9) [n=%d]"), c(nbrOfLoci(cnR), nbrOfLoci(cnSa), nbrOfLoci(cnSb))), bty="n")


# When target loci are identical to input ones, and unique,
# then the binned output should equal the input signals.
# Extract CNs with unique loci (most likely already the case)
x <- getPositions(cnR)
idxs <- which(!duplicated(x))
cnU <- cnR[idxs,]
xOut <- getPositions(cnU)
cnUs <- binnedSmoothing(cnU, xOut=xOut)
stopifnot(all(getPositions(cnUs) == getPositions(cnU)))
stopifnot(all(getSignals(cnUs) == getSignals(cnU)))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Binned smoothing (by count)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(cn, col="#999999", ylim=c(-3,3))
title(main="Binned smoothing (by count)")

cnSa <- binnedSmoothing(cn, by=3, byCount=TRUE)
lines(cnSa, col="blue")
points(cnSa, col="blue")

cnSb <- binnedSmoothing(cn, by=9, byCount=TRUE)
lines(cnSb, col="red")
points(cnSb, col="red")

legend("topright", pch=19, col=c("#999999", "blue", "red"), sprintf(c("raw [n=%d]", "BinO(w=3) [n=%d]", "BinO(w=9) [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cnSa), nbrOfLoci(cnSb))), bty="n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Kernel smoothing (default is Gaussian)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(cn, col="#999999", ylim=c(-3,3))
title(main="Kernel smoothing w/ Gaussian kernel")

cnSa <- kernelSmoothing(cn, h=2)
points(cnSa, col="blue")

cnSb <- kernelSmoothing(cn, h=5)
points(cnSb, col="red")

legend("topright", pch=19, col=c("#999999", "blue", "red"), sprintf(c("raw [n=%d]", "N(.,sd=2) [n=%d]", "N(.,sd=5) [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cnSa), nbrOfLoci(cnSb))), bty="n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Kernel smoothing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(cn, col="#999999", ylim=c(-3,3))
title(main="Kernel smoothing w/ uniform kernel")

xOut <- xSeq(cn, by=10)
cnSa <- kernelSmoothing(cn, xOut=xOut, kernel="uniform", h=2)
lines(cnSa, col="blue")
points(cnSa, col="blue")

cnSb <- kernelSmoothing(cn, xOut=xOut, kernel="uniform", h=5)
lines(cnSb, col="red")
points(cnSb, col="red")

legend("topright", pch=19, col=c("#999999", "blue", "red"), sprintf(c("raw [n=%d]", "U(w=2) [n=%d]", "U(w=5) [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cnSa), nbrOfLoci(cnSb))), bty="n")
