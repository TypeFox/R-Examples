.optimalWeights.simpleMV <- function(mu, sigma, constraints=NULL, tol = 1e-6)
{
	if(is.null(constraints))        
	{    
		numAssets <- length(mu)
		Amat <- rbind(rep(1, numAssets), diag(length(mu)))
		constraints <- list("Amat" = t(Amat), "bvec" = NULL, "meq" = 1)
		constraints$bvec <- c(1, rep(0, length(mu)))
		
	}
	stopifnot(class(constraints) == "list")
	stopifnot(all(c("Amat", "bvec", "meq") %in% names(constraints)))
	
	
	wts <- solve.QP(sigma, mu, constraints$Amat, constraints$bvec, constraints$meq)
#    else
#        wts <- solve.QP(sigma, mu, constraints$Amat, meq = constraints$meq)
	wts$solution[abs(wts$solution) < tol] <- 0
	names(wts$solution) <- names(mu)
	wts$solution
}



## A utility function that calculates "optimal" portfolios with respect to a prior and (Black-Litterman) posterior distribution, 
## and then returns the weights and optionally plots them with barplots. The optimizer is provided by the user, but there is a "toy" 
## Markowitz optimizer that is used by default. 

optimalPortfolios <- function
(                                         
	result, optimizer = .optimalWeights.simpleMV, ...,	doPlot = TRUE, 	beside = TRUE  
) 
{
	BARWIDTH <- 1
	
	.assertClass(result, "BLResult")
	optimizer <- match.fun(optimizer)
	
	# calculate the optimal prior and posterior weigths
	priorPortfolioWeights <- optimizer(result@priorMean, result@priorCovar, ...)
	postPortfolioWeights <- optimizer(result@posteriorMean, result@posteriorCovar, ...)
	if(doPlot)
	{        
		if(beside) {                                              
			plotData <- .removeZeroColumns(rbind(prior = priorPortfolioWeights, posterior = postPortfolioWeights))
			barplot(plotData, beside = TRUE,col = c("lightblue", "cyan"), border = "blue",
					legend.text = c("Prior", "Posterior"), horiz = FALSE, ylab = "Weights", main = "Optimal weights")
		}
		else
		{
			plotData <-  postPortfolioWeights -  priorPortfolioWeights
			plotData <- plotData[plotData != 0]
			barplot(plotData, col = c("lightblue"), ylab = "Difference", border = "blue", main = "Differences in weights", horiz = FALSE)
		}
		
	}
	return(list(priorPfolioWeights = priorPortfolioWeights, postPfolioWeights = postPortfolioWeights ))    
}

## A utility function that calculates "optimal" portfolios with respect to a prior and (Black-Litterman) 
## posterior distribution using the functionality of the Rmetrics fPortfolio package, and then returns the weights

optimalPortfolios.fPort <- function(result, spec = NULL, constraints = "LongOnly", optimizer = "minriskPortfolio", 
			inputData = NULL, numSimulations = BLCOPOptions("numSimulations"))
{
	stop("Not implemented for this class of result")
}

setGeneric("optimalPortfolios.fPort")

optimalPortfolios.fPort.BL <- function(result, spec = NULL ,constraints = "LongOnly", optimizer = "minriskPortfolio", 
		inputData = NULL, numSimulations = BLCOPOptions("numSimulations"))
{
#	if(!require("fPortfolio")) stop("The package fPortfolio is required to execute this function, but you do not have it installed.")
	assets <- assetSet(result@views)
	# create a "dummy" series that will only be used because the "optimizer" function requires it (but the mean and
	# covariance will not be calculated using it)
	dmySeries <- as.timeSeries(matrix(0, nrow = 1, ncol = length(assets), dimnames = list(NULL, assets)))
	numAssets <- length(assets)
	if(is.null(spec))
	{
		spec <- portfolioSpec()
		# # setType(spec) <- "MV"
		# setWeights(spec) <- rep(1 / numAssets, times = numAssets)
		 #setSolver(spec) <- "solveRquadprog"
	}
  
  # calculate prior and posterior mean and covariance. These are then stored in the scope of the
  # package environment and then accessed via a wrapper function when called by the optimiser.
  # This replaces the original code which stored value in the global environment and was therefore
  # in breach of the CRAN policies.
  
  .BLEnv$prior <- list(mu = result@priorMean, Sigma = result@priorCovar)
  
	# posterior mean and covariance estimates come from the BL calculations
	postMeanCov <- posteriorMeanCov(result)
	.BLEnv$post <- list(mu = postMeanCov$mean, Sigma = postMeanCov$covariance)
    
 	priorSpec <- spec
 	posteriorSpec <- spec
  
  setEstimator(priorSpec) <- "getPriorEstim"  
	setEstimator(posteriorSpec) <- "getPosteriorEstim"
  
  optimizer <- match.fun(optimizer)

  # calculate optimal portfolios
	priorOptimPortfolio <- optimizer(dmySeries, priorSpec, constraints)	
	posteriorOptimPortfolio <- optimizer(dmySeries, posteriorSpec, constraints)
  
	.BLEnv$prior <- NULL
	.BLEnv$post <-  NULL
  
	x <- list(priorOptimPortfolio = priorOptimPortfolio, posteriorOptimPortfolio = posteriorOptimPortfolio)
	class(x) <- "BLOptimPortfolios"
	x
}

## A wrapper function which returns estimates of the prior mean and covariance calculated 
## stored in the package environment. It is not intended to be directly run by the user but 
## needs to be exported in order to be run by the third party optimizer.

getPriorEstim <- function(x, spec=NULL, ...)
{
  return(.BLEnv$prior)
}


## A wrapper function which returns estimates of the posterior mean and covariance calculated 
## stored in the package environment. It is not intended to be directly run by the user but 
## needs to be exported in order to be run by the third party optimizer.

getPosteriorEstim <- function(x, spec=NULL, ...)
{
  return(.BLEnv$post)
}

setMethod("optimalPortfolios.fPort", signature(result = "BLResult"), optimalPortfolios.fPort.BL )

# plot methods, not yet exposed

plot.BLOptimPortfolios <- function(x,...)
{
	plotData <-   getWeights(x[[2]]@portfolio) - getWeights(x[[1]]@portfolio)
	
	plotData <- plotData[plotData != 0]
	
	barplot(plotData, col = c("lightblue"), ylab = "Difference", border = "blue", main = "Differences in weights", horiz = FALSE)	
}

## A utility function that calculates "optimal" portfolios with respect to a prior and (COP) 
## posterior distribution using the functionality of the Rmetrics fPortfolio package, and then returns the weights

optimalPortfolios.fPort.COP <- function(result, spec = NULL, constraints = "LongOnly", optimizer = "minriskPortfolio",
			inputData = NULL, numSimulations = BLCOPOptions("numSimulations"))
{

#	if(!require("fPortfolio")) stop("The package fPortfolio is required to execute this function, but you do not have it installed.")
	if(is.null(inputData))
	{
		# no input time series provided, so simulate the asset returns
		inputData <- sampleFrom(result@marketDist, n = numSimulations)
		colnames(inputData) <- assetSet(result@views)
		inputData <- as.timeSeries(inputData)
	}
	numAssets <- length(assetSet(result@views))
	# missing portfolio spec, create mean-CVaR portfolio
	if(is.null(spec))
	{
		spec <- portfolioSpec()
		setType(spec) <- "CVAR"
		setWeights(spec) <- rep(1 / numAssets, times = numAssets)
		setSolver(spec) <- "solveRglpk.CVAR"
	}
	# the "output" series data for the CVaR optimization will be taken from the simulations
	outData <- tail(posteriorSimulations(result), numSimulations)
	colnames(outData) <- assetSet(result@views)
	# kill row names to prevent errors when coercing to timeSeries object
	rownames(outData) <- NULL
	outData <- as.timeSeries(outData)
	optimizer <- match.fun(optimizer)
	# calculate prior and posterior optimal portfolios
	priorOptimPortfolio <- optimizer(inputData,spec, constraints)
	posteriorOptimPortfolio <- optimizer(outData, spec, constraints)
	x <- list(priorOptimPortfolio = priorOptimPortfolio, posteriorOptimPortfolio = posteriorOptimPortfolio)
	class(x) <- "COPOptimPortfolios"
	x
}

plot.COPOptimPortfolios <- function(x,...)
{
	plotData <-   getWeights(x[[2]]@portfolio) - getWeights(x[[1]]@portfolio)
	
	plotData <- plotData[plotData != 0]
	
	barplot(plotData, col = c("lightblue"), ylab = "Difference", border = "blue", main = "Differences in weights", horiz = FALSE)	
}


setMethod("optimalPortfolios.fPort", signature(result = "COPResult"), optimalPortfolios.fPort.COP )
