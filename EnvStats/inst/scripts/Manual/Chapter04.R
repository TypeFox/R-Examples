#####
# File:     Chapter04.R
#
# Purpose:  Reproduce Examples in Chapter 4 of the book:
#
#           Millard, SP. (2013).  EnvStats: An R Package for 
#             Environmental Statistics.  Springer-Verlag.
#
# Author:   Steven P. Millard
#
# Last
# Updated:  2013/02/03
#####


#######################################################################################

library(EnvStats)

############################################################

# 4.1 Introduction 
#-----------------

# Figure 4.1
#-----------

# Figure 4.1a
#------------
windows()
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
pdfPlot(dist = "beta", param.list = list(shape1=2, shape2=4),
  curve.fill.col = "cyan")
pdfPlot(dist = "beta", param.list = list(shape1=1, shape2=1, ncp=1),
  curve.fill.col = "cyan")
pdfPlot(dist = "binom", param.list = list(size=10, prob=0.5),
  hist.col = "cyan")
pdfPlot(dist = "cauchy", param.list = list(location=0, scale=1),
  left.tail.cutoff = 0.01, right.tail.cutoff = 0.01, 
  curve.fill.col = "cyan")
pdfPlot(dist = "chi", param.list = list(df=4),
  curve.fill.col = "cyan")
pdfPlot(dist = "chisq", param.list = list(df=4),
  curve.fill.col = "cyan")


# Figure 4.1b
#------------
windows()
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
pdfPlot(dist = "chisq", param.list = list(df=5, ncp=1),
  curve.fill.col = "cyan")
set.seed(21)
epdfPlot(rgamma(100, shape=4, scale=5), curve.fill = TRUE, curve.fill.col = "cyan", 
  xlab = "Observations", 
  main = "Empirical Density Based On 100\nGamma(shape=4, scale=5) Random Numbers",
  cex.main = 1)
pdfPlot(dist = "exp", param.list = list(rate=2),
  curve.fill.col = "cyan")
pdfPlot(dist = "evd", param.list = list(location=0, scale=1),
  curve.fill.col = "cyan")
pdfPlot(dist = "gevd", param.list = list(location=0, scale=1, shape = 0.5),
  curve.fill.col = "cyan", cex.main = 1)
pdfPlot(dist = "f", param.list = list(df1=5, df2=10),
  curve.fill.col = "cyan")


# Figure 4.1c
#------------
windows()
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
pdfPlot(dist = "f", param.list = list(df1=5, df2=10, ncp=1),
  curve.fill.col = "cyan")
pdfPlot(dist = "gamma", param.list = list(shape=2, scale=1),
  curve.fill.col = "cyan")
pdfPlot(dist = "gammaAlt", param.list = list(mean=10, cv=0.5),
  curve.fill.col = "cyan")
pdfPlot(dist = "geom", param.list = list(prob=0.5),
  hist.col = "cyan")
pdfPlot(dist = "hyper", param.list = list(m=20, n=15, k=7),
  hist.col = "cyan")
pdfPlot(dist = "logis", param.list = list(location=0, scale=1),
  curve.fill.col = "cyan")


# Figure 4.1d
#------------
windows()
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
pdfPlot(dist = "lnorm", param.list = list(meanlog=0, sdlog=1),
  curve.fill.col = "cyan")
pdfPlot(dist = "lnormAlt", param.list = list(mean=10, cv=0.5),
  curve.fill.col = "cyan")
pdfPlot(dist = "lnormMix", param.list = 
  list(meanlog1=0, sdlog1=1, meanlog2=3, sdlog2=0.5, p.mix=0.5),
  right.tail.cutoff = 0.02,
  curve.fill.col = "cyan", cex.main = 1, 
  main = paste("Lognormal Mixture Density", "(meanlog1=0, sdlog1=1,",
  "meanlog2=3, sdlog2=0.5, p.mix=0.5)", sep="\n"))
pdfPlot(dist = "lnormMixAlt", param.list = 
  list(mean1=5, cv1=1, mean2=20, cv2=0.5, p.mix=0.5),
  right.tail.cutoff = 0.01,
  curve.fill.col = "cyan", cex.main=1,
  main = paste("Lognormal Mixture Density", "(mean1=5, cv1=1,",
  "mean2=20, cv2=0.5, p.mix=0.5)", sep="\n"))
pdfPlot(dist = "lnorm3", param.list = list(meanlog=0, sdlog=1, threshold=5),
  right.tail.cutoff = 0.01, curve.fill.col = "cyan", cex.main = 1)
pdfPlot(dist = "lnormTrunc", param.list = 
  list(meanlog=0, sdlog=1, min=0, max=2),
  curve.fill.col = "cyan", cex.main = 1)


# Figure 4.1e
#------------
windows()
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
pdfPlot(dist = "lnormTruncAlt", param.list = 
  list(mean=2, cv=1, min=0, max=3),
  curve.fill.col = "cyan")
pdfPlot(dist = "nbinom", param.list = list(size=4, prob=0.5),
  hist.col = "cyan")
pdfPlot(dist = "norm", param.list = list(mean=0, sd=1),
  curve.fill.col = "cyan")
pdfPlot(dist = "normMix", param.list = 
  list(mean1=0, sd1=1, mean2=4, sd2=2, p.mix=0.5),
  curve.fill.col = "cyan", cex.main=1,
  main = paste("Normal Mixture Density", "(mean1=0, sd1=1,",
  "mean2=4, sd2=2, p.mix=0.5)", sep="\n"))
pdfPlot(dist = "normTrunc", param.list = 
  list(mean=10, sd=2, min=8, max=13),
  curve.fill.col = "cyan", cex.main = 1)
pdfPlot(dist = "pareto", param.list = list(location=1, shape=2),
  curve.fill.col = "cyan", right.tail.cutoff = 0.01)


# Figure 4.1f
#------------
windows()
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
pdfPlot(dist = "pois", param.list = list(lambda=5),
  hist.col = "cyan")
pdfPlot(dist = "t", param.list = list(df=5),
  curve.fill.col = "cyan")
pdfPlot(dist = "t", param.list = list(df=5, ncp=1),
  curve.fill.col = "cyan")
pdfPlot(dist = "tri", param.list = list(min=0, max=1, mode=0.7),
  curve.fill.col = "cyan")
pdfPlot(dist = "unif", param.list = list(min=0, max=1),
  curve.fill.col = "cyan")
pdfPlot(dist = "weibul", param.list = list(shape=2, scale=1),
  curve.fill.col = "cyan")


# Figure 4.1g
#------------
windows()
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), mgp = c(1.5, 0.5, 0))
pdfPlot(dist = "wilcox", param.list = list(m=4, n=3),
  hist.col = "cyan")
pdfPlot(dist = "zmlnorm", param.list = list(meanlog=0, sdlog=1, p.zero=0.5),
  right.tail.cutoff = 0.01, curve.fill.col = "cyan", cex.main = 1)
pdfPlot(dist = "zmlnormAlt", param.list = list(mean=2, cv=1, p.zero=0.4),
  right.tail.cutoff = 0.01, curve.fill.col = "cyan", cex.main = 1)
pdfPlot(dist = "zmnorm", param.list = list(mean=5, sd=1, p.zero=0.3),
  curve.fill.col = "cyan")
frame()
frame()



# 4.2.1 Probability Density Function for Lognormal Distribution
#--------------------------------------------------------------

# Figure 4.2
#-----------

windows()
with(EPA.94b.tccb.df, 
  hist(TcCB[Area == "Reference"], freq = FALSE, 
  xlim = c(0, 2), xlab = "TcCB (ppb)", col = "cyan",
  main = "Density Histogram of Reference Area TcCB Data"))


# Figure 4.3
#-----------
windows()
pdfPlot(distribution = "lnormAlt", param.list = list(mean = 0.6, cv = 0.5), 
  curve.fill.col = "cyan")


round(dlnormAlt(seq(0, 2, by = 0.5), mean = 0.6, cv = 0.5), 3)

# Figure 4.4
#-----------
windows()
pdfPlot(distribution = "gammaAlt", param.list = list(mean = 0.6, cv = 0.5))

round(dgammaAlt(seq(0, 2, by = 0.5), mean = 0.6, cv = 0.5), 3)


# 4.3.1 Cumulative Distribution Function for Lognormal Distribution
#------------------------------------------------------------------

# Figure 4.5
#-----------
windows()
cdfPlot(distribution = "lnormAlt", param.list = list(mean = 0.6, cv = 0.5))

round(plnormAlt(seq(0, 2, by = 0.5), mean = 0.6, cv = 0.5), 2)


# 4.4.1 Quantiles for Lognormal Distribution
#-------------------------------------------
qlnormAlt(c(0.5, 0.95), mean = 0.6, cv = 0.5)



# 4.5.1 Generating Random Numbers from a Univariate Distribution
#---------------------------------------------------------------
set.seed(23)
rlnormAlt(5, mean = 0.6, cv = 0.5)


# 4.5.2 Generating Multivariate Normal Random Numbers
#----------------------------------------------------
library(MASS)
set.seed(47)
sd.vec <- c(1, 3)
cor.mat <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
cov.mat <- diag(sd.vec) %*% cor.mat %*% diag(sd.vec)
mvrnorm(n = 3, mu = c(5, 10), Sigma = cov.mat)

rm(sd.vec, cor.mat, cov.mat)


# 4.5.3 Generating Multivariate Observations Based on Rank Correlations
#----------------------------------------------------------------------
simulateMvMatrix(n = 3, 
  distributions = c(X1 = "norm", X2 = "lnormAlt"),
  param.list = list(X1 = list(mean = 5, sd = 1),
                    X2 = list(mean = 10, cv = 2)),
  cor.mat = matrix(c(1, 0.5, 0.5, 1), ncol=2), seed = 105)


