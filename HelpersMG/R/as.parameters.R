#' @title Extract parameters from mcmcComposite object
#' @author Marc Girondot
#' @return A vector with parameters at maximum likelihood or index position
#' @param x A mcmcComposite obtained as a result of \code{MHalgoGen()} function
#' @param index At which iteration the parameters must be taken
#' @param chain The number of the chain inwhich to get parameters
#' @description Take a mcmcComposite object and create a vector object with parameter value at specified iteration.\cr
#' If \code{index="best"}, the function will return the parameters for the highest likelihood. It also indicates at which iteration the maximum lihelihood has been observed.\cr
#' If \code{index="last"}, the fuction will return the parameters for the last likelihood.\cr
#' \code{index} can also be a numeric value.
#' @family mcmcComposite functions
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' require(coda)
#' x <- rnorm(30, 10, 2)
#' dnormx <- function(x, par) return(-sum(dnorm(x, mean=par['mean'], sd=par['sd'], log=TRUE)))
#' parameters_mcmc <- data.frame(Density=c('dnorm', 'dlnorm'), 
#' Prior1=c(10, 0.5), Prior2=c(2, 0.5), SDProp=c(1, 1), 
#' Min=c(-3, 0), Max=c(100, 10), Init=c(10, 2), stringsAsFactors = FALSE, 
#' row.names=c('mean', 'sd'))
#' mcmc_run <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=1)
#' plot(mcmc_run, xlim=c(0, 20))
#' plot(mcmc_run, xlim=c(0, 10), parameters="sd")
#' mcmcforcoda <- as.mcmc(mcmc_run)
#' #' heidel.diag(mcmcforcoda)
#' raftery.diag(mcmcforcoda)
#' autocorr.diag(mcmcforcoda)
#' acf(mcmcforcoda[[1]][,"mean"], lag.max=20, bty="n", las=1)
#' acf(mcmcforcoda[[1]][,"sd"], lag.max=20, bty="n", las=1)
#' batchSE(mcmcforcoda, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmcforcoda)$statistics[,"Time-series SE"]
#' summary(mcmc_run)
#' as.parameters(mcmc_run)
#' lastp <- as.parameters(mcmc_run, index="last")
#' parameters_mcmc[,"Init"] <- lastp
#' # The n.adapt set to 1 is used to not record the first set of parameters
#' # then it is not duplicated (as it is also the last one for 
#' # the object mcmc_run)
#' mcmc_run2 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=1, thin=1, trace=1)
#' mcmc_run3 <- merge(mcmc_run, mcmc_run2)
#' ####### no adaptation, n.adapt must be 0
#' parameters_mcmc[,"Init"] <- c(mean(x), sd(x))
#' mcmc_run3 <- MHalgoGen(n.iter=1000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=0, thin=1, trace=1)
#' }
#' @export


as.parameters <-
function(x, index="best", chain=1) {

	L <- x$resultLnL[[chain]]
	p <- x$resultMCMC[,][[chain]]
	
	tab <- 0
	pos <- NULL
	
	if (is.numeric(index)) {
	  pos <- index
	} else {
	
if (index=="best") {
  pos <- which.max(L)
  print(paste("The best likelihood has been observed at iteration", pos-tab))
}

if (index=="last") {
  pos <- length(L)
}
	}
	
	if (is.null(pos)) {
	  warning("index is not recognized")
	  return()
	}
	
	pml <- as.numeric(p[pos,])
	names(pml) <- names(p[pos,])
	return(pml)

}
