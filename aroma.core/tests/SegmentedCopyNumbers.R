library("aroma.core")
set.seed(0xBEEF)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# True CN states
stateFcn <- function(x, ...) {
  states <- integer(length(x))
  states[100 <=x & x <= 200] <- -1L
  states[320 <=x & x <= 400] <- +1L
  states
}

# Number of loci
J <- 500

x <- 1:J
y <- rnorm(J, sd=1/2)

# Reshuffle
o <- sample(J)
x <- x[o]
y <- y[o]

for (state in c(-1,+1)) {
  idxs <- (stateFcn(x) == state)
  y[idxs] <- y[idxs] + state
}

cn <- SegmentedCopyNumbers(y, x, states=stateFcn)
print(cn)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sorting
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnO <- sort(cn)
print(cnO)

# Sanity check
oinv <- order(o)
stopifnot(all.equal(cn[oinv,], cnO))
o <- order(cn$x)
stopifnot(identical(o, oinv))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subsetting
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnS <- extractSubsetByState(cn, states=c(0,+1L))

# Sanity check pre and post sorting
cnOS <- extractSubsetByState(cnO, states=c(0,+1L))
stopifnot(all.equal(cnOS, sort(cnS)))

plot(cn, ylim=c(-4,4))
title("Copy numbers annotated by state (and subset by state)")
print(cnS)
points(cnS, pch=21, cex=1.2, lwd=2, col="purple")
legend("topright", pch=c(19, 21), col=c("#999999", "purple"), sprintf(c("raw [n=%d]", "CN in {0,1} [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cnS))), bty="n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Binned smoothing stratified by state
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnSa <- binnedSmoothingByState(cn, by=3, verbose=-1)
cnSb <- binnedSmoothingByState(cn, by=9, verbose=-1)

# Sanity check pre and post sorting
cnOSb <- binnedSmoothingByState(cnO, by=9)
stopifnot(all.equal(cnOSb, sort(cnSb)))

plot(cn, col="#999999", ylim=c(-3,3))
title(main="Binned smoothing stratified by state")
lines(cnSa, col="blue")
points(cnSa, col="blue")
lines(cnSb, col="red")
points(cnSb, col="red")
legend("topright", pch=19, col=c("#999999", "blue", "red"), sprintf(c("raw [n=%d]", "Bin(w=3) [n=%d]", "Bin(w=9) [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cnSa), nbrOfLoci(cnSb))), bty="n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Kernel smoothing stratified by state
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnSa <- kernelSmoothingByState(cn, h=2)
cnSb <- kernelSmoothingByState(cn, h=5)

# Sanity check pre and post sorting
cnOSb <- kernelSmoothingByState(cnO, h=5)
##stopifnot(all.equal(cnOSb, sort(cnSb)))

plot(cn, col="#999999", ylim=c(-3,3))
title(main="Kernel smoothing stratified by state w/ Gaussian kernel")
points(cnSa, col="blue")
points(cnSb, col="red")
legend("topright", pch=19, col=c("#999999", "blue", "red"), sprintf(c("raw [n=%d]", "N(.,sd=2) [n=%d]", "N(.,sd=5) [n=%d]"), c(nbrOfLoci(cn), nbrOfLoci(cnSa), nbrOfLoci(cnSb))), bty="n")
