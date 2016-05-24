library("aroma.core")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data from multiple platforms
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Piecewise-constant copy-number state function
cnState <- function(x, ...) {
  n <- length(x)
  mu <- double(n)
  mu[20e6 <= x & x <= 30e6] <- +1
  mu[65e6 <= x & x <= 80e6] <- -1
  mu
} # cnState()

xMax <- 100e6

Js <- c(200, 400, 100)
bs <- c(1, 1.4, 0.5)
as <- c(0, +0.5, -0.5)
sds <- c(0.5, 0.3, 0.8)

cnList <- list()
for (kk in seq(along=Js)) {
  J <- Js[kk]
  a <- as[kk]
  b <- bs[kk]
  sd <- sds[kk]
  x <- sort(runif(J, max=xMax))
  mu <- a + b * cnState(x)
  eps <- rnorm(J, sd=sd)
  y <- mu + eps
  cn <- RawCopyNumbers(y, x)
  cnList[[kk]] <- cn
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Merge platform data (record their sources in 'id')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cn <- Reduce(append, cnList)
plot(cn, ylim=c(-3,3), col=cn$id)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
legend <- c()

if (require("DNAcopy")) {
  fit <- segmentByCBS(cn)
  cnr <- extractCopyNumberRegions(fit)
  print(cnr)
  drawLevels(cnr, col="white", lwd=6)
  drawLevels(cnr, col="red", lwd=3)
  legend <- c(legend, red="CBS")
}


if (require("mpcbs")) {
  fit <- segmentByMPCBS(cn)
  cnr <- extractCopyNumberRegions(fit)
  print(cnr)
  drawLevels(cnr, col="white", lwd=6)
  drawLevels(cnr, col="blue", lwd=3)
  legend <- c(legend, blue="MPCBS")
}


if (length(legend) > 0) {
  legend("topleft", pch=19, col=names(legend), legend, bty="n", horiz=TRUE)
}
