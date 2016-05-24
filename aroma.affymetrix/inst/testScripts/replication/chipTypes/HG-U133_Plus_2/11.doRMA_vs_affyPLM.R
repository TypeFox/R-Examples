###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the RMA
# chip-effect estimates as estimated by affyPLM.
#
# Author: Mark Robinson and Henrik Bengtsson
# Created: 2007-06-20
# Last modified: 2012-09-01
###########################################################################
library("aroma.affymetrix")
library("affyPLM")          # fitPLM()
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

# Avoid being masked by oligo::fitPLM() and affy::plotDensity()
fitPLM <- affyPLM::fitPLM
plotDensity <- aroma.light::plotDensity

# Detach 'oligoClasses' in case it is loaded.  If not, there an error
# related to probeNames() will be thrown.
tryCatch(detach("package:oligoClasses"), error=function(ex) {})


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
verbose && enter(verbose, "RMA by aroma.affymetrix")

res <- doRMA(csR, drop=FALSE, verbose=verbose)
print(res)

# Extract chip effects on the log2 scale
ces <- res$ces
theta <- extractMatrix(ces)
rownames(theta) <- getUnitNames(cdf)
theta <- log2(theta)

verbose && exit(verbose)


# ------------------------
# RMA estimates by affyPLM
# ------------------------
verbose && enter(verbose, "RMA by affyPLM")
verbose && print(verbose, sessionInfo())

raw <- ReadAffy(filenames=getPathnames(csR))
verbose && print(verbose, raw)

fit <- fitPLM(raw, verbos=9)
theta0 <- coefs(fit)

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
toPNG(getFullName(csR), tags=c("doRMA_vs_affyPLM"), width=800, {
  par(mar=c(5,5,4,2)+0.1, cex.main=2, cex.lab=2, cex.axis=1.5)

  layout(matrix(1:16, ncol=4, byrow=TRUE))

  xlab <- expression(log[2](theta[affyPLM]))
  ylab <- expression(log[2](theta[aroma.affymetrix]))
  for (kk in seq_len(ncol(theta))) {
    main <- colnames(theta)[kk]
    plot(theta0[,kk], theta[,kk], pch=".", xlab=xlab, ylab=ylab, main=main)
    abline(0,1, col="blue")
    stext(side=3, pos=0, line=-1.1, cex=1.2, substitute(rho==x, list(x=rho[kk])))
  }

  xlab <- expression(log[2](theta[aroma.affymetrix]/theta[affyPLM]))
  plotDensity(e, xlab=xlab)
})

# (b) Assert correlations
print(rho)
print(range(rho))
stopifnot(all(rho > 0.99995))

# (c) Assert differences
stopifnot(mean(as.vector(e^2)) < 1e-3)
stopifnot(sd(as.vector(e^2)) < 1e-3)
stopifnot(quantile(abs(e), 0.99) < 0.05)
stopifnot(max(abs(e)) < 0.100)


verbose && print(verbose, sessionInfo())


###########################################################################
# HISTORY:
# 2012-09-01 [HB]
# o Now using a GEO data set.
# 2012-08-30 [HB]
# o Updated to utilize toPNG().
# 2010-02-05 [HB]
# o Harmonized with the corresponding gcRMA test script.
# 2010-01-02 [HB]
# o BUG FIX: If loaded, detaching oligoClasses, because otherwise it will
#   cause name conflict with affy::probeNames().
# 2008-07-17 [HB]
# o Added some more quantile-based assertions too.
# o Had to lower the similarity threshold from 1e-4 to 1e-3. I don't know
#   why this is, but the differences are still very small.
###########################################################################
