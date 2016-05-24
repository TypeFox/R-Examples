#' phenology_MHmcmc runs the Metropolis-Hastings algorithm for data (Bayesian MCMC)
#' @title Run the Metropolis-Hastings algorithm for data
#' @author Marc Girondot
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param n.iter Number of iterations for each step
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param result An object obtained after a SearchR fit
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param trace True or False, shows progress
#' @param filename If intermediate is not NULL, save intermediate result in this file
#' @param intermediate Period for saving intermediate result, NULL for no save
#' @param previous Previous result to be continued. Can be the filename in which intermediate results are saved.
#' @description Run the Metropolis-Hastings algorithm for data.\cr
#' The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
#' I recommend thin=1 because the method to estimate SE uses resampling.\cr
#' As initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' The parameters intermediate and filename are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file of name filename.\cr
#' The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
#' These options are used to prevent the consequences of computer crash or if the run is very very long and computer processes at time limited.\cr
#' @examples 
#' \dontrun{
#' library(phenology)
#' data(Gratiot)
#' # Generate a formatted list named data_Gratiot 
#' data_Gratiot <- add_phenology(Gratiot, name="Complete", 
#'     reference=as.Date("2001-01-01"), format="%d/%m/%Y")
#' # Generate initial points for the optimisation
#' parg <- par_init(data_Gratiot, parametersfixed=NULL)
#' # Run the optimisation
#' result_Gratiot <- fit_phenology(data=data_Gratiot, 
#' 		parametersfit=parg, parametersfixed=NULL, trace=1)
#' # Generate set of priors for Bayesian analysis
#' pmcmc <- phenology_MHmcmc_p(result_Gratiot, accept = TRUE)
#' result_Gratiot_mcmc <- phenology_MHmcmc(result = result_Gratiot, n.iter = 10000, 
#' parametersMCMC = pmcmc, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' # Get standard error of parameters
#' summary(result_Gratiot_mcmc)
#' # Make diagnostics of the mcmc results using coda package
#' mcmc <- as.mcmc(result_Gratiot_mcmc)
#' require(coda)
#' heidel.diag(mcmc)
#' raftery.diag(mcmc)
#' autocorr.diag(mcmc)
#' acf(mcmc[[1]][,"LengthB"], lag.max=200, bty="n", las=1)
#' acf(mcmc[[1]][,"Max_Gratiot"], lag.max=50, bty="n", las=1)
#' batchSE(mcmc, batchSize=100)
#' # The batch standard error procedure is usually thought to 
#' # be not as accurate as the time series methods used in summary
#' summary(mcmc)$statistics[,"Time-series SE"]
#' plot(result_Gratiot_mcmc, parameters=3, las=1, xlim=c(-10, 300))
#' }
#' @export


phenology_MHmcmc<-function(result=stop("An output from fit_phenology() must be provided"), n.iter=10000, 
parametersMCMC=stop("A model generated with phenology_MHmcmc_p() must be provided"), n.chains = 4, 
n.adapt = 0, thin=1, trace=FALSE, intermediate=NULL, filename="intermediate.Rdata", previous=NULL) {
  
  # result <- result_Gratiot; n.iter <- 10000;parametersMCMC <- pmcmc; n.chains = 1; n.adapt = 0; thin = 1; trace = FALSE
  # intermediate=NULL; filename="intermediate.Rdata"; previous=NULL
  
  if (is.character(previous)) {
    itr <- NULL
    load(previous)
    previous <- itr
    rm(itr)
    print("Continue previous mcmc run")
  }

  if (class(result)!="phenology") {
    stop("An output from fit_phenology() must be provided")
  }
  

# likelihood <- .phenology_fonctionMCMC

pt=list(data=result$data, fixed=result$parametersfixed, incertitude=result$method_incertitude, zerocounts=result$zero_counts, infinite=result$infinite, out=TRUE)

print(parametersMCMC)

out <- MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, n.chains = n.chains, n.adapt = n.adapt, thin=thin, trace=trace, pt=pt, likelihood=getFromNamespace(".Lnegbin", ns="phenology"))

fin <- try(summary(out), silent=TRUE)

if (class(fin)=="try-error") {
  lp <- rep(NA, nrow(out$parametersMCMC$parameters))
  names(lp) <- rownames(out$parametersMCMC$parameters)
  out <- c(out, SD=list(lp))
} else {
  out <- c(out, SD=list(fin$statistics[,"SD"]))
}

class(out) <- "mcmcComposite"

return(out)

}
