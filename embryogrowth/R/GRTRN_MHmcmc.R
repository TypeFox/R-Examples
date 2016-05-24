#' GRTRN_MHmcmc runs the Metropolis-Hastings algorithm for data (Bayesian MCMC)
#' @title Metropolis-Hastings algorithm for Embryo Growth Rate Thermal Reaction Norm
#' @author Marc Girondot
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param n.iter Number of iterations for each step
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param result An object obtained after a SearchR fit
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param trace True or False, shows progress
#' @param batchSize Number of observations to include in each batch fo SE estimation
#' @param parallel If true, try to use several cores using parallel computing
#' @param filename If intermediate is not NULL, save intermediate result in this file
#' @param intermediate Period for saving intermediate result, NULL for no save
#' @param previous Previous result to be continued. Can be the filename in which intermediate results are saved.
#' @description Run the Metropolis-Hastings algorithm for data.\cr
#' The number of iterations is \code{n.iter+n.adapt+1} because the initial likelihood is also displayed.\cr
#' I recommend that thin=1 because the method to estimate SE uses resampling.\cr
#' If initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' To get the SE of the point estimates from \code{result_mcmc <- GRTRN_MHmcmc(result=try)}, use:\cr
#' \code{result_mcmc$SD}\cr
#' \code{coda} package is necessary for this function.\cr
#' The parameters \code{intermediate} and \code{filename} are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file named \code{filename}.\cr
#' The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
#' These options are used to prevent the consequences of computer crash or if the run is very very long and processes with user limited time.\cr
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(nest)
#' formated <- FormatNests(nest)
#' # The initial parameters value can be:
#' # "T12H", "DHA",  "DHH", "Rho25"
#' # Or
#' # "T12L", "T12H", "DHA",  "DHH", "DHL", "Rho25"
#' x <- structure(c(118.768297442004, 475.750095909406, 306.243694918151, 
#' 116.055824800264), .Names = c("DHA", "DHH", "T12H", "Rho25"))
#' # pfixed <- c(K=82.33) or rK=82.33/39.33
#' pfixed <- c(rK=2.093313)
#' resultNest_4p <- searchR(parameters=x, fixed.parameters=pfixed,  
#' 	temperatures=formated, derivate=dydt.Gompertz, M0=1.7,  
#' 	test=c(Mean=39.33, SD=1.92))
#' data(resultNest_4p)
#' pMCMC <- TRN_MHmcmc_p(resultNest_4p, accept=TRUE)
#' # Take care, it can be very long; several days
#' result_mcmc_4p <- GRTRN_MHmcmc(result=resultNest_4p, 
#' 		parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' 		n.adapt = 0, thin=1, trace=TRUE)
#' data(result_mcmc_4p)
#' out <- as.mcmc(result_mcmc_4p)
#' # This out can be used with coda package
#' # Test for stationarity and length of chain
#' require(coda)
#' heidel.diag(out)
#' raftery.diag(out)
#' # plot() can use the direct output of GRTRN_MHmcmc() function.
#' plot(result_mcmc_4p, parameters=1, xlim=c(0,550))
#' plot(result_mcmc_4p, parameters=3, xlim=c(290,320))
#' # summary() permits to get rapidly the standard errors for parameters
#' # They are store in the result also.
#' se <- result_mcmc_4p$SD
#' # the confidence interval is better estimated by:
#' apply(out[[1]], 2, quantile, probs=c(0.025, 0.975))
#' # The use of the intermediate method is as followed;
#' # Here the total mcmc iteration is 10000, but every 1000, intermediate
#' # results are saved in file intermediate1000.Rdata:
#' result_mcmc_4p <- GRTRN_MHmcmc(result=resultNest_4p, 
#' parametersMCMC=pMCMC, n.iter=10000, n.chains = 1,  
#' n.adapt = 0, thin=1, trace=TRUE, 
#' intermediate=1000, filename="intermediate1000.Rdata")
#' # If run has been stopped for any reason, it can be resumed with:
#' result_mcmc_4p <- GRTRN_MHmcmc(previous="intermediate1000.Rdata")
#' }
#' @export

GRTRN_MHmcmc <- function(result=NULL, n.iter=10000, 
parametersMCMC=NULL, n.chains = 1, n.adapt = 0, 
thin=1, trace=FALSE, batchSize=sqrt(n.iter), parallel=TRUE, 
intermediate=NULL, filename="intermediate.Rdata", previous=NULL)
{

  # result=NULL; n.iter=10000; parametersMCMC=NULL; n.chains = 1; n.adapt = 0; thin=1; trace=FALSE; batchSize=sqrt(n.iter); parallel=TRUE; intermediate=NULL; filename="intermediate.Rdata"; previous=NULL
  # previous <- "/Users/marc/Dropbox/DropBoxPerso/_Package_HelpersMG/previous_result_mcmc_AtlMed_simplify_2.RData"
  # result=resultNest_4p; parametersMCMC=pMCMC; n.iter=100; n.chains = 1; n.adapt = 0; thin=1; trace=TRUE; intermediate = 10
  if (is.character(previous)) {
    itr <- NULL
    load(previous)
    previous <- itr
    rm(itr)
  }
  
  if (is.list(previous)) {
    print("Continue previous mcmc run")
    # 29/1/2014; Ajout de result$weight
    out <- MHalgoGen(previous=previous)
    
  } else {
    print(parametersMCMC)
    # 29/1/2014; Ajout de result$weight
    # n.iter=n.iter; parameters=parametersMCMC; n.chains = n.chains; n.adapt = n.adapt; thin=thin; trace=trace; data=list(data=result$data, derivate=result$derivate, test=result$test, M0=result$M0, fixed.parameters=result$fixed.parameters, weight=result$weight); likelihood=getFromNamespace(".fonctionMCMC", ns="embryogrowth"); intermediate=intermediate; filename=filename; previous=previous
    # 3/12/2015 j'avais un data=list(data=result$data, derivate=result$derivate, 
    # test=result$test, M0=result$M0, fixed.parameters=result$fixed.parameters, 
    # weight=result$weight)
    out <- MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, 
                                  n.chains = n.chains, n.adapt = n.adapt, thin=thin, trace=trace, 
                                  data=result$data, derivate=result$derivate, 
                                            test=result$test, M0=result$M0, fixed.parameters=result$fixed.parameters, 
                                            weight=result$weight, 
                                  likelihood=getFromNamespace(".fonctionMCMC", ns="embryogrowth"), 
                                  intermediate=intermediate, filename=filename, previous=previous)
    
  }


if (batchSize>=n.iter/2) {
  print("batchSize cannot be larger than half the number of iterations.")
  rese <- rep(NA, dim(parametersMCMC)[1])
  names(rese) <- rownames(parametersMCMC)
  out <- c(out, SE=list(rese))
} else {
  out <- c(out, BatchSE=list(coda::batchSE(out$resultMCMC, batchSize=batchSize)))
}

class(out) <- "mcmcComposite"

fin <- try(summary(out), silent=TRUE)

if (class(fin)=="try-error") {
  lp <- rep(NA, nrow(out$parametersMCMC$parameters))
  names(lp) <- rownames(out$parametersMCMC$parameters)
  out <- c(out, TimeSeriesSE=list(lp))
  out <- c(out, SD=list(lp))
} else {
  out <- c(out, TimeSeriesSE=list(fin$statistics[,4]))
  out <- c(out, SD=list(fin$statistics[,"SD"]))
}

class(out) <- "mcmcComposite"

return(out)
}
