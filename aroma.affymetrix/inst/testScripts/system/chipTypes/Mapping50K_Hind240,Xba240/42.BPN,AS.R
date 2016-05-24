library("aroma.affymetrix")
library("matrixStats"); # rowMedians()

# Avoid being masked by affy::plotDensity()
plotDensity <- aroma.light::plotDensity

log <- Arguments$getVerbose(-4, timestamp=TRUE)

dataSet <- "HapMap,CEU,testset"
chipType <- "Mapping50K_Hind240"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup of annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CDF
cdf <- AffymetrixCdfFile$byChipType(chipType)

# Assert existence of probe-sequence annotation files
acs <- AromaCellSequenceFile$byChipType(chipType)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf)
print(csR)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR)
print(acc)
csC <- process(acc, verbose=log)
print(csC)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Base-position normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csC, shift=+300)
print(bpn)

csN <- process(bpn, verbose=log)
print(csN)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allele-specific chip effect estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE)
print(plm)
fit(plm, verbose=log)
ces <- getChipEffectSet(plm)
print(ces)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces)
print(fln)
cesN <- process(fln, verbose=log)
print(cesN)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract (thetaA, thetaB) for copy-neutral chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- getCdf(cesN)
gi <- getGenomeInformation(cdf)
units <- getUnitsOnChromosomes(gi, 1:22)
theta <- extractTheta(cesN, units=units)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate copy numbers on the natural scale
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
thetaR <- rowMedians(theta[,1,]+theta[,2,], na.rm=TRUE)
C <- 2*(theta[,1,]+theta[,2,])/thetaR
CA <- 2*theta[,1,]/thetaR
CB <- 2*theta[,2,]/thetaR
freqB <- CB/C

mu <- colMedians(C, na.rm=TRUE)
muA <- colMedians(CA, na.rm=TRUE)
muB <- colMedians(CB, na.rm=TRUE)

Clim <- c(-0.2,3.2)
xlab <- "Freq B (thetaB/theta)"
Clab <- "Copy number"
CAlab <- "Copy number (Allele A)"
CBlab <- "Copy number (Allele B)"
col <- "green"

fig <- 1
if (!devIsOpen(fig <- fig + 1)) {
  devNew()
  layout(matrix(1:9, nrow=3, byrow=TRUE))
  par(mar=c(4,4,0.5,0.5)+0.1)
  centers <- matrix(c(0,2, 1,1, 2,0), nrow=3, ncol=2, byrow=TRUE)
  for (cc in 1:ncol(C)) {
    xx <- CA[,cc]
    yy <- CB[,cc]
    X <- cbind(xx,yy)
    ok <- (is.finite(X) & -1 < X & X < 30)
    ok <- ok[,1] & ok[,2]
    X <- X[ok,]
    smoothScatter(X, xlim=Clim, ylim=Clim, xlab=CAlab, ylab=CBlab)
    abline(h=0:3, lty=3, lwd=1)
    abline(v=0:3, lty=3, lwd=1)
    abline(a=2, b=-1, lty=1, lwd=1, col="white")
    abline(a=2, b=-1, lty=2, lwd=1)

    # Plot centers
    fit <- kmeans(X, centers=centers, algorithm="Lloyd", iter.max=30)
    print(fit$centers)
    lines(fit$centers[c(1,3),], col=col, lwd=1)
    points(fit$centers, pch=21, col="white", bg=col, lwd=1.5)
  }
  devDone()
}

if (!devIsOpen(fig <- fig + 1)) {
  devNew()
  layout(matrix(1:9, nrow=3, byrow=TRUE))
  par(mar=c(4,4,0.5,0.5)+0.1)
  centers <- matrix(c(0,2, 1/2,2, 1,2), nrow=3, ncol=2, byrow=TRUE)
  for (cc in 1:ncol(C)) {
    xx <- freqB[,cc]
    yy <- C[,cc]
    X <- cbind(xx,yy)
    ok <- (is.finite(X) & -1 < X & X < 30)
    ok <- ok[,1] & ok[,2]
    X <- X[ok,]
    smoothScatter(X, xlim=c(0,1), ylim=Clim, xlab=xlab, ylab=Clab)
    abline(h=0:3, lty=3, lwd=1)
    abline(v=0:2/2, lty=3, lwd=1)

    # Plot centers
    fit <- kmeans(X, centers=centers, algorithm="Lloyd", iter.max=30)
    print(fit$centers)
    lines(fit$centers[c(1,3),], col=col, lwd=1)
    points(fit$centers, pch=21, col="white", bg=col, lwd=1.5)
  }
  plotDensity(freqB, lwd=2, xlim=c(0,1), ylim=c(0,2.2))
  abline(v=c(0,1/2,1), lty=3, lwd=1)
  abline(h=c(0,1,2), lty=3, lwd=1)
  devDone()
}
