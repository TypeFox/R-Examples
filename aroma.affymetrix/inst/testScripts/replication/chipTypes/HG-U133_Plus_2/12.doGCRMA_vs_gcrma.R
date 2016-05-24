###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the gcRMA
# chip-effect estimates as estimated by gcrma.
#
# Author: Mark Robinson and Henrik Bengtsson
# Created: 2009-05-17
# Last modified: 2012-09-01
###########################################################################
library("aroma.affymetrix")
library("gcrma")  # gcrma()
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

# Avoid being masked by affy::plotDensity()
plotDensity <- aroma.light::plotDensity


# ----------------------------------
# Dataset
# ----------------------------------
dataSet <- "GSE9890"
chipType <- "HG-U133_Plus_2"

cdf <- AffymetrixCdfFile$byChipType(chipType)
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf)
print(csR)


# ----------------------------------
# RMA estimates by aroma.affymetrix
# ----------------------------------
verbose && enter(verbose, "gcRMA by aroma.affymetrix")

res <- doGCRMA(csR, drop=FALSE, verbose=verbose)
print(res)

# Extract chip effects on the log2 scale
ces <- res$ces
theta <- extractMatrix(ces)
rownames(theta) <- getUnitNames(cdf)
theta <- log2(theta)

verbose && exit(verbose)


# ------------------------
# gcRMA estimates by gcrma
# ------------------------
verbose && enter(verbose, "gcRMA by gcrma")
verbose && print(verbose, sessionInfo())

raw <- ReadAffy(filenames=getPathnames(csR))
verbose && print(verbose, raw)

es <- gcrma(raw, verbose=TRUE)
verbose && print(verbose, es)

theta0 <- exprs(es)
verbose && exit(verbose)


# --------------------------------
# Compare the two implementations
# --------------------------------
# Reorder the aroma.affymetrix estimates
o <- match(rownames(theta0), rownames(theta))
theta <- theta[o,]

# Calculate statistics
rho <- diag(cor(theta, theta0))
print(rho)
print(range(rho))
e <- (theta - theta0)
print(summary(e))

# (a) Visual comparison
toPNG(getFullName(csR), tags=c("doGCRMA_vs_gcrma"), width=800, {
  par(mar=c(5,5,4,2)+0.1, cex.main=2, cex.lab=2, cex.axis=1.5)

  layout(matrix(1:16, ncol=4, byrow=TRUE))

  xlab <- expression(log[2](theta[gcrma]))
  ylab <- expression(log[2](theta[aroma.affymetrix]))
  for (kk in seq_len(ncol(theta))) {
    main <- colnames(theta)[kk]
    plot(theta0[,kk], theta[,kk], pch=".", xlab=xlab, ylab=ylab, main=main)
    abline(0,1, col="blue")
    stext(side=3, pos=0, line=-1.1, cex=1.2, substitute(rho==x, list(x=rho[kk])))
  }

  xlab <- expression(log[2](theta[aroma.affymetrix]/theta[gcrma]))
  plotDensity(e, xlab=xlab)
})

# (b) Assert correlations
stopifnot(all(rho > 0.9999))

# (c) Assert differences
stopifnot(mean(as.vector(e^2)) < 0.001)
stopifnot(sd(as.vector(e^2)) < 0.003)
stopifnot(quantile(abs(e), 0.99) < 0.10)
stopifnot(max(abs(e)) < 0.30)


verbose && print(verbose, sessionInfo())

###########################################################################
# HISTORY:
# 2012-09-01 [HB]
# o Updated to use a GEO data set.
# 2012-08-30 [HB]
# o Updated to utilize toPNG().
# 2010-02-05 [HB]
# o Harmonized with the corresponding RMA test script.
# 2009-05-17 [HB]
# o Adopted from Mark Robinson's test script.
###########################################################################
