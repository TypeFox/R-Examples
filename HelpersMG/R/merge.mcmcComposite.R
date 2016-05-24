#' merge.mcmcComposite Merge two mcmcComposite results
#' @title Merge two mcmcComposite results
#' @author Marc Girondot
#' @return A mcmcComposite result
#' @param x A mcmcComposite obtained as a result of \code{MHalgoGen()} function
#' @param y A mcmcComposite obtained as a result of \code{MHalgoGen()} function
#' @param ... not used
#' @family mcmcComposite functions
#' @description Merge two mcmcComposite results and produced a new one mcmcComposite object.\cr
#' Note that the initial value 
#' for the second run must use the last value of the first one as shown in example.
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
#' @method merge mcmcComposite
#' @export

merge.mcmcComposite <- function(x, y, ...) {
  
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("coda package is necessary for this function")
  }
  
  
  mcmcComposite1 <- x
  mcmcComposite2 <- y
  
  if (mcmcComposite1$parametersMCMC$n.chains != mcmcComposite2$parametersMCMC$n.chains) {
    print("MCMC results must have the same number of chains")
    return()
  }
  
  mcmcComposite <- mcmcComposite1
  for (chain in 1:mcmcComposite1$parametersMCMC$n.chains) {
  mcmcComposite$resultMCMC[[chain]] <- mcmc(rbind(mcmcComposite1$resultMCMC[[chain]], 
                                             mcmcComposite2$resultMCMC[[chain]]), 1, 
                                            mcmcComposite1$parametersMCMC$n.iter+
                                              mcmcComposite2$parametersMCMC$n.iter-
                                              mcmcComposite1$parametersMCMC$n.adapt-
                                              mcmcComposite2$parametersMCMC$n.adapt, 1)
  }
  
  mcmcComposite$resultLnL[[1]] <- c(mcmcComposite1$resultLnL[[1]], 
                               mcmcComposite2$resultLnL[[1]])
  
  mcmcComposite$parametersMCMC$n.iter <- mcmcComposite1$parametersMCMC$n.iter+
    mcmcComposite2$parametersMCMC$n.iter
  
  e <- mcmcComposite$resultMCMC
  
  mcmcComposite$BatchSE <- batchSE(e)
  mcmcComposite$TimeSeriesSE <- summary(e)$statistics[,"Time-series SE"]
  return(mcmcComposite)
}
