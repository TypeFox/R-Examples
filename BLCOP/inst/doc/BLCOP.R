### R code from vignette source 'BLCOP.Rnw'

###################################################
### code chunk number 1: BLCOP.Rnw:76-79
###################################################
require(fPortfolio)
require(BLCOP)
require(mnormt)


###################################################
### code chunk number 2: BLCOP.Rnw:81-85
###################################################
pickMatrix <- matrix(c(1/2, -1, 1/2, rep(0, 3)), nrow = 1, ncol = 6 )
views <- BLViews(P = pickMatrix, q = 0.06,confidences =  100, 
assetNames = colnames(monthlyReturns))
views


###################################################
### code chunk number 3: BLCOP.Rnw:90-92
###################################################
priorMeans <- rep(0, 6)
priorVarcov <- cov.mve(monthlyReturns)$cov


###################################################
### code chunk number 4: BLCOP.Rnw:100-102
###################################################
marketPosterior <- posteriorEst(views = views, sigma = priorVarcov, 
 mu = priorMeans, tau = 1/2)


###################################################
### code chunk number 5: BLCOP.Rnw:107-111
###################################################
finViews <- matrix(ncol = 4, nrow = 1, dimnames = list(NULL, c("C","JPM","BAC","MS")))
finViews[,1:4] <- rep(1/4,4)
views <- addBLViews(finViews, 0.15, 90, views)
views


###################################################
### code chunk number 6: BLCOP.Rnw:118-120
###################################################
marketPosterior <- BLPosterior(as.matrix(monthlyReturns), views, tau = 1/2, 
	marketIndex = as.matrix(sp500Returns),riskFree = as.matrix(US13wTB))


###################################################
### code chunk number 7: BLCOP.Rnw:136-140
###################################################
optPorts <- optimalPortfolios.fPort(marketPosterior, optimizer = "tangencyPortfolio")
par(mfcol = c(2, 1))
weightsPie(optPorts$priorOptimPortfolio)
weightsPie(optPorts$posteriorOptimPortfolio)


###################################################
### code chunk number 8: BLCOP.Rnw:146-149
###################################################
optPorts2 <- optimalPortfolios.fPort(marketPosterior, 
		constraints = "minW[1:6]=0.1", optimizer = "minriskPortfolio")
optPorts2


###################################################
### code chunk number 9: BLCOP.Rnw:155-156
###################################################
densityPlots(marketPosterior, assetsSel = "JPM")


###################################################
### code chunk number 10: BLCOP.Rnw:207-214
###################################################
dispersion <- c(.376,.253,.360,.333,.360,.600,.397,.396,.578,.775) / 1000
sigma <- BLCOP:::.symmetricMatrix(dispersion, dim = 4)
caps <- rep(1/4, 4)
mu <- 2.5 * sigma %*% caps
dim(mu) <- NULL
marketDistribution <- mvdistribution("mt", mean = mu, S = sigma, df = 5 )
class(marketDistribution)


###################################################
### code chunk number 11: BLCOP.Rnw:222-228
###################################################
pick <- matrix(0, ncol = 4, nrow = 1, 
	dimnames = list(NULL, c("SP", "FTSE", "CAC", "DAX")))
pick[1,"DAX"] <- 1
viewDist <- list(distribution("unif", min = -0.02, max = 0))
views <- COPViews(pick, viewDist = viewDist, confidences = 0.2, 
	assetNames = c("SP", "FTSE", "CAC", "DAX"))


###################################################
### code chunk number 12: BLCOP.Rnw:233-238
###################################################
newPick <- matrix(0, 1, 2)
dimnames(newPick) <- list(NULL, c("SP", "FTSE"))
newPick[1,] <- c(1, -1) # add a relative view
views <- addCOPViews(newPick, 
	list(distribution("norm", mean = 0.05, sd = 0.02)), 0.5, views)


###################################################
### code chunk number 13: BLCOP.Rnw:243-245
###################################################
marketPosterior <- COPPosterior(marketDistribution, views, numSimulations = 50000)
densityPlots(marketPosterior, assetsSel = 4)


