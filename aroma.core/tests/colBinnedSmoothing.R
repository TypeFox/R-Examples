library("aroma.core")

# Number of tracks
I <- 4

# Number of data points per track
J <- 100

# Simulate data with a gain in track 2 and 3
x <- 1:J
Y <- matrix(rnorm(I*J, sd=1/2), ncol=I)
Y[30:50,2:3] <- Y[30:50,2:3] + 3

# Uniformly distributed equal-sized bins
Ys3 <- colBinnedSmoothing(Y, x=x, from=2, by=3)
Ys5 <- colBinnedSmoothing(Y, x=x, from=3, by=5)

# Custom bins
xOutRange <- t(matrix(c(
  1, 11,
 11, 31,
 31, 41,
 41, 51,
 51, 81,
 81, 91,
 91,101
), nrow=2))
YsC <- colBinnedSmoothing(Y, x=x, xOutRange=xOutRange)

# Custom bins specified by center locations with
# maximized width relative to the neighboring bins.
xOut <- c(6, 21, 36, 46, 66, 86, 96)
YsD <- colBinnedSmoothing(Y, x=x, xOut=xOut)

xlim <- range(x)
ylim <- c(-3,5)
layout(matrix(1:I, ncol=1))
par(mar=c(3,3,1,1)+0.1, pch=19)
for (ii in 1:I) {
  plot(NA, xlim=xlim, ylim=ylim)
  points(x, Y[,ii], col="#999999")

  xOut <- attr(Ys3, "xOut")
  lines(xOut, Ys3[,ii], col=2)
  points(xOut, Ys3[,ii], col=2)

  xOut <- attr(Ys5, "xOut")
  lines(xOut, Ys5[,ii], col=3)
  points(xOut, Ys5[,ii], col=3)

  xOut <- attr(YsC, "xOut")
  lines(xOut, YsC[,ii], col=4)
  points(xOut, YsC[,ii], col=4, pch=15)

  xOut <- attr(YsD, "xOut")
  lines(xOut, YsD[,ii], col=5)
  points(xOut, YsD[,ii], col=5, pch=15)

  if (ii == 1) {
    legend("topright", pch=c(19,19,15,15), col=c(2,3,4,5),
           c("by=3", "by=5", "Custom #1", "Custom #2"), horiz=TRUE, bty="n")
  }
}


# Sanity checks
xOut <- x
YsT <- colBinnedSmoothing(Y, x=x, xOut=xOut)
stopifnot(all(YsT == Y))
stopifnot(all(attr(YsT, "counts") == 1))

xOut <- attr(YsD, "xOut")
YsE <- colBinnedSmoothing(YsD, x=xOut, xOut=xOut)
stopifnot(all(YsE == YsD))
stopifnot(all(attr(YsE, "xOutRange") == attr(YsD, "xOutRange")))
stopifnot(all(attr(YsE, "counts") == 1))

# Scramble ordering of loci
idxs <- sample(x)
x2 <- x[idxs]
Y2 <- Y[idxs,,drop=FALSE]
Y2s <- colBinnedSmoothing(Y2, x=x2, xOut=x2)
stopifnot(all(attr(Y2s, "xOut") == x2))
stopifnot(all(attr(Y2s, "counts") == 1))
stopifnot(all(Y2s == Y2))

xOut <- x[seq(from=2, to=J, by=3)]
YsT <- colBinnedSmoothing(Y, x=x, xOut=xOut)
stopifnot(all(YsT == Ys3))
stopifnot(all(attr(YsT, "counts") == 3))

xOut <- x[seq(from=3, to=J, by=5)]
YsT <- colBinnedSmoothing(Y, x=x, xOut=xOut)
stopifnot(all(YsT == Ys5))
stopifnot(all(attr(YsT, "counts") == 5))

