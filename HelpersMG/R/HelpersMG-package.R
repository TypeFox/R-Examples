#' Helpers functions for several packages
#'
#' \tabular{ll}{
#'  Package: \tab HelpersMG\cr
#'  Type: \tab Package\cr
#'  Version: \tab 1.5 build 51\cr
#'  Date: \tab 2016-05-01\cr
#'  License: \tab GPL (>= 2)\cr
#'  LazyLoad: \tab yes\cr
#'  }
#' @title Tools for Earth Meteorological Analysis
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType package
#' @name HelpersMG-package
#' @description Contains functions useful for managing \cr
#' 'NetCDF' files (see http://en.wikipedia.org/wiki/NetCDF), \cr
#' get tide levels on any point of the globe, \cr
#' get moon phase and time for sun rise and fall, \cr
#' analyse and reconstruct daily time series of temperature \cr
#' with irregular sinusoidal pattern, \cr
#' show scales and wind rose in plot with change of color of text, \cr 
#' Metropolis-Hastings algorithm for Bayesian MCMC analysis, \cr
#' plot graphs or boxplot with error bars, \cr
#' search files in disk by there names or their content, \cr
#' read the contents of all files from a folder at one time, \cr
#' calculate IC50 for ecotoxicological studies, \cr
#' calculate the probability mass function of the sum of negative binomial \cr
#' distributions.\cr
#' The lastest version of this package can always been installed using:\cr
#' install.packages("http://www.ese.u-psud.fr/epc/conservation/CRAN/HelpersMG.tar.gz", 
#' repos=NULL, type="source")
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' print('----------------------------------')
#' print('Examples for mcmcComposite objects')
#' print('----------------------------------')
#' require(coda)
#' x <- rnorm(30, 10, 2)
#' dnormx <- function(x, par) return(-sum(dnorm(x, mean=par['mean'], sd=par['sd'], log=TRUE)))
#' parameters_mcmc <- data.frame(Density=c('dnorm', 'dlnorm'), 
#' Prior1=c(10, 0.5), Prior2=c(2, 0.5), SDProp=c(0.35, 0.2), 
#' Min=c(-3, 0), Max=c(100, 10), Init=c(10, 2), stringsAsFactors = FALSE, 
#' row.names=c('mean', 'sd'))
#' mcmc_run <- MHalgoGen(n.iter=100000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=100, thin=1, trace=1)
#' plot(mcmc_run, xlim=c(0, 20))
#' plot(mcmc_run, xlim=c(0, 10), parameters="sd")
#' mcmcforcoda <- as.mcmc(mcmc_run)
#' # Optimal rejection rate should be 0.234
#' rejectionRate(mcmcforcoda)
#' heidel.diag(mcmcforcoda)
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
#' mcmc_run2 <- MHalgoGen(n.iter=10000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=1, thin=1, trace=1)
#' mcmc_run3 <- merge(mcmc_run, mcmc_run2)
#' ####### no adaptation, n.adapt must be 0
#' parameters_mcmc[,"Init"] <- c(mean(x), sd(x))
#' mcmc_run3 <- MHalgoGen(n.iter=10000, parameters=parameters_mcmc, data=x, 
#' likelihood=dnormx, n.chains=1, n.adapt=0, thin=1, trace=1)
#' print('------------------------------------------')
#' print('Examples for Daily patterns of temperature')
#' print('------------------------------------------')
#' # Generate a timeserie of time
#' time.obs <- NULL
#' for (i in 0:9) time.obs <- c(time.obs, c(0, 6, 12, 18)+i*24)
#' # For these time, generate a timeseries of temperatures
#' temp.obs <- rep(NA, length(time.obs))
#' temp.obs[3+(0:9)*4] <- rnorm(10, 25, 3)
#' temp.obs[1+(0:9)*4] <- rnorm(10, 10, 3)
#' for (i in 1:(length(time.obs)-1)) 
#'   if (is.na(temp.obs[i])) 
#'   temp.obs[i] <- mean(c(temp.obs[i-1], temp.obs[i+1]))
#'   if (is.na(temp.obs[length(time.obs)])) 
#'   temp.obs[length(time.obs)] <- temp.obs[length(time.obs)-1]/2
#' observed <- data.frame(time=time.obs, temperature=temp.obs)
#' # Search for the minimum and maximum values
#' r <- minmax.periodic(time.minmax.daily=c(Min=2, Max=15), 
#' observed=observed, period=24)
#' 
#' # Estimate all the temperatures for these values
#' t <- temperature.periodic(minmax=r)
#' 
#' plot_errbar(x=t[,"time"], y=t[,"temperature"],
#' errbar.y=ifelse(is.na(t[,"sd"]), 0, 2*t[,"sd"]),
#' type="l", las=1, bty="n", errbar.y.polygon = TRUE, 
#' xlab="hours", ylab="Temperatures", ylim=c(0, 35), 
#' errbar.y.polygon.list = list(col="grey"))
#' 
#' plot_add(x=t[,"time"], y=t[,"temperature"], type="l")
#' }

NULL

