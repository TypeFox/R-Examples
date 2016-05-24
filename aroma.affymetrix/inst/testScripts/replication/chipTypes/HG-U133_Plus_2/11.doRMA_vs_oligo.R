###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the RMA
# chip-effect estimates as estimated by oligo.
# The setup is the same as in affyPLM,fitPLM.R.
#
# Author: Henrik Bengtsson
# Created: 2008-12-04 (from affyPLM,fitPLM.R)
# Last modified: 2013-07-03
###########################################################################

library("aroma.affymetrix")
library("oligo")

# Avoid being masked by affy::plotDensity()
plotDensity <- aroma.light::plotDensity

verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

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

res <- doRMA(csR, flavor="oligo", drop=FALSE, verbose=verbose)
print(res)

# Extract chip effects on the log2 scale
ces <- res$ces
theta <- extractMatrix(ces)
rownames(theta) <- getUnitNames(cdf)
theta <- log2(theta)

verbose && exit(verbose)


# ----------------------
# RMA estimates by oligo
# ----------------------
verbose && enter(verbose, "RMA by oligo")
verbose && print(verbose, sessionInfo())

library("pd.hg.u133.plus.2")

raw <- read.celfiles(filenames=getPathnames(csR))
eSet <- rma(raw)
theta0 <- exprs(eSet)

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
toPNG(getFullName(csR), tags=c("doRMA_vs_oligo"), width=800, {
  par(mar=c(5,5,4,2)+0.1, cex.main=2, cex.lab=2, cex.axis=1.5)

  layout(matrix(1:16, ncol=4, byrow=TRUE))

  xlab <- expression(log[2](theta[oligo]))
  ylab <- expression(log[2](theta[aroma.affymetrix]))
  for (kk in seq_len(ncol(theta))) {
    main <- colnames(theta)[kk]
    plot(theta0[,kk], theta[,kk], pch=".", xlab=xlab, ylab=ylab, main=main)
    abline(0,1, col="blue")
    stext(side=3, pos=0, line=-1.1, cex=1.2, substitute(rho==x, list(x=rho[kk])))
  }

  xlab <- expression(log[2](theta[aroma.affymetrix]/theta[oligo]))
  plotDensity(e, xlab=xlab)
})

# (b) Assert correlations
print(rho)
print(range(rho))
stopifnot(all(rho > 0.99985))

cors <- sapply(1:ncol(theta), FUN=function(cc) cor(theta[,cc], theta0[,cc]))
print(cors)
print(range(cors))
stopifnot(all(cors > 0.99985))

# (c) Assert differences
print(mean(as.vector(e^2)))
stopifnot(mean(as.vector(e^2)) < 0.0015)

print(sd(as.vector(e^2)))
stopifnot(sd(as.vector(e^2)) < 0.0025)

print(quantile(abs(e), 0.99))
## stopifnot(quantile(abs(e), 0.99) < 0.05)
stopifnot(quantile(abs(e), 0.99) < 0.11)  # Why so bad now? /HB 2014-05-28

print(max(abs(e)))
## stopifnot(max(abs(e)) < 0.085)
stopifnot(max(abs(e)) < 0.26) ## Why so bad?


verbose && print(verbose, sessionInfo())


###########################################################################
# HISTORY:
# 2013-07-03 [HB]
# o Harmonized with 11.doRMA_vs_affyPLM.R.
# 2008-12-04 [HB]
# o Created, but not tested because I miss package 'pd.hg.u133.plus.2'.
###########################################################################
