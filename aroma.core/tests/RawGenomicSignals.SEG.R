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
x <- sort(runif(length(y), max=length(y))) * 1e5;
w <- runif(J)
w[320:400] <- 0.001


cn <- RawCopyNumbers(y, x)
print(cn)

plot(cn, ylim=c(-3,3), col="#aaaaaa", xlab="Position (Mb)")

cnS <- binnedSmoothing(cn, by=500e3)
print(cnS)
lines(cnS, col="black", lwd=3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
legend <- c()

if (require("DNAcopy")) {
  fit <- segmentByCBS(cn)
  cnr <- extractCopyNumberRegions(fit)
  print(cnr)
  drawLevels(cnr, col="red", lwd=3)
  legend <- c(legend, red="CBS")
}

if (require("GLAD")) {
  fit <- segmentByGLAD(cn)
  cnr <- extractCopyNumberRegions(fit)
  print(cnr)
  drawLevels(cnr, col="blue", lwd=3)
  legend <- c(legend, blue="GLAD")
}

if (require("HaarSeg")) {
  fit <- segmentByHaarSeg(cn)
  cnr <- extractCopyNumberRegions(fit)
  print(cnr)
  drawLevels(cnr, col="orange", lwd=3)
  legend <- c(legend, orange="HaarSeg")
}

if (require("mpcbs")) {
  fit <- segmentByMPCBS(cn)
  cnr <- extractCopyNumberRegions(fit)
  print(cnr)
  drawLevels(cnr, col="white", lwd=6)
  drawLevels(cnr, col="purple", lwd=3)
  legend <- c(legend, purple="MPCBS")
}
 

if (length(legend) > 0) {
  legend("topleft", pch=19, col=names(legend), legend, bty="n", horiz=TRUE)
}
