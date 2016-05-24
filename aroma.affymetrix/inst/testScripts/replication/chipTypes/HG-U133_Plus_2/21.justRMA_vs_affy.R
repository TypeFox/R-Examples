###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the justRMA()
# chip-effect estimates as estimated by the affy package.
#
# Author: Henrik Bengtsson
# Created: 2014-04-28
###########################################################################
library("aroma.affymetrix")
library("affy")
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
# affy
# ----------------------------------
verbose && enter(verbose, "justRMA() by aroma.affymetrix")
eset0 <- justRMA(filenames=getPathnames(csR))
print(eset0)
verbose && exit(verbose)


# ----------------------------------
# aroma.affymetrix
# ----------------------------------
verbose && enter(verbose, "justRMA() by aroma.affymetrix")
eset <- justRMA(csR, verbose=verbose)
print(eset)
verbose && exit(verbose)


# --------------------------------
# Compare the two implementations
# --------------------------------
# Sample names
sampleNames0 <- sampleNames(eset0)
sampleNames <- sampleNames(eset)
stopifnot(identical(sampleNames, tools::file_path_sans_ext(sampleNames0)))

# Feature names
stopifnot(identical(featureNames(eset), featureNames(eset0)))

# Annotation
stopifnot(identical(annotation(eset), annotation(eset0)))

# Expression estimates
theta0 <- exprs(eset0)
theta <- exprs(eset)

# Calculate statistics
rho <- diag(cor(theta, theta0))
print(rho)
print(range(rho))
e <- (theta - theta0)
print(summary(e))

# (a) Visual comparison
toPNG(getFullName(csR), tags=c("justRMA", "aroma.affymetrix_vs_affy"), width=800, {
  par(mar=c(5,5,4,2)+0.1, cex.main=2, cex.lab=2, cex.axis=1.5)

  layout(matrix(1:16, ncol=4, byrow=TRUE))

  xlab <- expression(log[2](theta[affy]))
  ylab <- expression(log[2](theta[aroma.affymetrix]))
  for (kk in seq_len(ncol(theta))) {
    main <- colnames(theta)[kk]
    plot(theta0[,kk], theta[,kk], pch=".", xlab=xlab, ylab=ylab, main=main)
    abline(0,1, col="blue")
    stext(side=3, pos=0, line=-1.1, cex=1.2, substitute(rho==x, list(x=rho[kk])))
  }

  xlab <- expression(log[2](theta[aroma.affymetrix]/theta[affy]))
  plotDensity(e, xlab=xlab)
})

# (b) Assert correlations
print(rho)
print(range(rho))
stopifnot(all(rho > 0.9998))

# (c) Assert differences
stopifnot(mean(as.vector(e^2)) < 0.002)
stopifnot(sd(as.vector(e^2)) < 0.003)
stopifnot(quantile(abs(e), 0.99) < 0.11)
stopifnot(max(abs(e)) < 0.35)


verbose && print(verbose, sessionInfo())


###########################################################################
# HISTORY:
# 2014-04-28 [HB]
# o Created from '11.doRMA_vs_affyPLM.R'.
###########################################################################
