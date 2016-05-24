# $LastChangedDate: 2010-02-28 13:45:31 +0000 (Sun, 28 Feb 2010) $
# $Rev: 4767 $
# Author: Francisco
###############################################################################
BLExample <- function()
{
	entries <- c(0.001005,0.001328,-0.000579,-0.000675,0.000121,0.000128,-0.000445,-0.000437 ,
			0.001328,0.007277,-0.001307,-0.000610,-0.002237,-0.000989,0.001442,-0.001535 ,
			-0.000579,-0.001307,0.059852,0.027588,0.063497,0.023036,0.032967,0.048039 ,
			-0.000675,-0.000610,0.027588,0.029609,0.026572,0.021465,0.020697,0.029854 ,
			0.000121,-0.002237,0.063497,0.026572,0.102488,0.042744,0.039943,0.065994 ,
			0.000128,-0.000989,0.023036,0.021465,0.042744,0.032056,0.019881,0.032235 ,
			-0.000445,0.001442,0.032967,0.020697,0.039943,0.019881,0.028355,0.035064 ,
			-0.000437,-0.001535,0.048039,0.029854,0.065994,0.032235,0.035064,0.079958 )
	
	myVarcov2 <- matrix(entries, ncol = 8, nrow = 8)
	mu <- c(0.08, 0.67,6.41, 4.08, 7.43, 3.70, 4.80, 6.60) / 100
	pick <- matrix(0, ncol = 8, nrow = 3, dimnames = list(NULL, letters[1:8]))
	pick[1,7] <- 1
	pick[2,1] <- -1; pick[2,2] <- 1
	pick[3, 3:6] <- c(0.9, -0.9, .1, -.1)
	confidences <- 1 / c(0.00709, 0.000141, 0.000866)
	myViews <- BLViews(pick, c(0.0525, 0.0025, 0.02), confidences, letters[1:8])
	myPosterior <- posteriorEst(myViews, tau = 0.025, mu = mu, myVarcov2 )
	list(prior = myViews, posterior = myPosterior)
}

# check the basic "toy" portfolio optimizer and the portfolio optimizer that uses fPortfolio
test.optimalPortfolios.BL <- function()
{
    BLEx <- get(load( file.path(BLCOPOptions("unitTestPath"), "BLExample.RData") ))
	#BLEx <- get(load(BLExample.RData))
	#BLEx <- BLCOP:::BLExample()
	myPosterior <- BLEx$posterior
	res <- optimalPortfolios(myPosterior, doPlot = TRUE)
	expected <- c(0.00000000,0.00000000,0.38204176, 0.00000000, 0.08198505,0.00000000, 0.36548138,0.17049181)
	checkEquals(res$priorPfolioWeights, expected, checkNames = FALSE, tolerance = 1e-06)
	expected <- c(0.00000000, 0.00000000, 0.32444891, 0.08071719, 0.09377903, 0.00000000 ,0.32478427 ,0.17627060)
	checkEquals(res$postPfolioWeights, expected, checkNames = FALSE, tolerance = 1e-06)
	
	# now try the fPortfolio optimizer
	
	#optimalPortfolios.fPort.BL <- function(result, spec,constraints = "LongOnly", optimizer = "minriskPortfolio", 
	#		inputData = NULL, numSimulations = NA)
	if(!require("fPortfolio"))
	{
		warning("The fPortfolio package is required to run these tests, but you don't have it installed")
		return()
	}
	res2 <- optimalPortfolios.fPort(myPosterior, spec = portfolioSpec(), optimizer = "tangencyPortfolio")
	
	# there should be two portfolios, each of class
	checkEquals(c(class(res2[[1]]), class(res2[[2]])), c("fPORTFOLIO", "fPORTFOLIO"), check.names = FALSE)
	
	# posterior weights should be similar to those in Idzorek's paper
	checkEquals(round(getWeights(res2$"posteriorOptimPortfolio"), 2), c(0.28, 0.17, 0.1, 0.14, 0.01, 0.01, 0.25, 0.04), check.names=FALSE)
	# try another optimizer
	
	res3 <- optimalPortfolios.fPort(myPosterior, spec = portfolioSpec(), optimizer = "minriskPortfolio")
	
	checkEquals(round(sapply(res3,  getWeights), 2), 
        structure(c(0.94, 0, 0, 0.04, 0, 0, 0.02, 0, 0.94, 0, 0, 0.04, 
        0, 0, 0.02, 0), .Dim = c(8L, 2L), .Dimnames = list(c("a", "b", 
        "c", "d", "e", "f", "g", "h"), c("priorOptimPortfolio", "posteriorOptimPortfolio"
        ))) ,msg = " |minimum risk portfolios as expected")
	
}

# tests optimalPortfolios.fPort for COPResults objects

test.optimalPortfolios.COP <- function()
{

	if(!require("fPortfolio", quiet = TRUE))
	{
		warning("This test relies on the fPortfolio package which is not available \n")
		return()
	}
	
	COPEx <- get(load( file.path(BLCOPOptions("unitTestPath"), "copexample.RData") ))
	# Check optimization with COP
	myPosterior <- COPEx$posterior
	
	res4 <- optimalPortfolios.fPort(myPosterior, spec = NULL, optimizer = "minriskPortfolio", inputData = NULL, 
			numSimulations  = 100	) 
	
	checkEqualsNumeric(round(getWeights(res4$posteriorOptimPortfolio),3),  c(0.535, 0.465, 0, 0) )
	
	# second example, using input data
	COPEx2 <- get(load( file.path(BLCOPOptions("unitTestPath"), "copexample2.RData") ))
	
	spec <- portfolioSpec()
	setType(spec) <- "CVaR"
	setWeights(spec) <- rep(1 / 6, times = 6)
	setSolver(spec) <- "solveRglpk.CVAR"
	setTargetReturn(spec) <- 0.005
	res5 <- optimalPortfolios.fPort( COPEx2, spec = spec, inputData = as.timeSeries(monthlyReturns), numSimulations  = nrow(monthlyReturns))

    checkEqualsNumeric(round(getWeights(res5$priorOptimPortfolio),3), c(0.071, 0, 0.011, 0.522, 0, 0.396))
                   
    checkEqualsNumeric(round(getWeights(res5$posteriorOptimPortfolio),3), c(0.513, 0, 0, 0, 0, 0.487))
}
