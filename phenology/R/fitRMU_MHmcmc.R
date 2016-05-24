#' fitRMU_MHmcmc runs the Metropolis-Hastings algorithm for RMU.data (Bayesian MCMC)
#' @title Run the Metropolis-Hastings algorithm for RMU.data
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
#' @family Fill gaps in RMU
#' @description Run the Metropolis-Hastings algorithm for RMU.data.\cr
#' The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
#' I recommend thin=1 because the method to estimate SE uses resampling.\cr
#' As initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' The parameters intermediate and filename are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file of name filename.\cr
#' The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
#' These options are used to prevent the consequences of computer crash or if the run is very very long and computer processes at time limited.\cr
#' @examples 
#' \dontrun{
#' library("phenology")
#' RMU.names.AtlanticW <- data.frame(mean=c("Yalimapo.French.Guiana", 
#'                                          "Galibi.Suriname", 
#'                                          "Irakumpapy.French.Guiana"), 
#'                                  se=c("se_Yalimapo.French.Guiana", 
#'                                       "se_Galibi.Suriname", 
#'                                       "se_Irakumpapy.French.Guiana"))
#' data.AtlanticW <- data.frame(Year=c(1990:2000), 
#'       Yalimapo.French.Guiana=c(2076, 2765, 2890, 2678, NA, 
#'                                6542, 5678, 1243, NA, 1566, 1566),
#'       se_Yalimapo.French.Guiana=c(123.2, 27.7, 62.5, 126, NA, 
#'                                  230, 129, 167, NA, 145, 20),
#'       Galibi.Suriname=c(276, 275, 290, NA, 267, 
#'                        542, 678, NA, 243, 156, 123),
#'       se_Galibi.Suriname=c(22.3, 34.2, 23.2, NA, 23.2, 
#'                            4.3, 2.3, NA, 10.3, 10.1, 8.9),
#'       Irakumpapy.French.Guiana=c(1076, 1765, 1390, 1678, NA, 
#'                                3542, 2678, 243, NA, 566, 566),
#'       se_Irakumpapy.French.Guiana=c(23.2, 29.7, 22.5, 226, NA, 
#'                                  130, 29, 67, NA, 15, 20))
#'                            
#' cst <- fitRMU(RMU.data=data.AtlanticW, RMU.names=RMU.names.AtlanticW, 
#'                colname.year="Year", model.trend="Constant", 
#'                model.SD="Zero")
#' pMCMC <- fitRMU_MHmcmc_p(result=cst, accept=TRUE)
#' fitRMU_MCMC <- fitRMU_MHmcmc(result = cst, n.iter = 10000, 
#' parametersMCMC = pMCMC, n.chains = 1, n.adapt = 0, thin = 1, trace = FALSE)
#' }
#' @export


fitRMU_MHmcmc <- function(result=stop("An output from fitRMU_MHmcmc() must be provided"), n.iter=10000, 
parametersMCMC=stop("A model generated with fitRMU_MHmcmc_p() must be provided"), n.chains = 4, 
n.adapt = 0, thin=1, trace=FALSE, intermediate=NULL, filename="intermediate.Rdata", previous=NULL) {
  
  if (is.character(previous)) {
    itr <- NULL
    load(previous)
    previous <- itr
    rm(itr)
    print("Continue previous mcmc run")
  }

  if (class(result)!="fitRMU") {
    warning("An output from fitRMU() must be provided")
    return()
  }

  fun <- getFromNamespace(".LikelihoodRMU", ns="phenology")

  print(parametersMCMC)
  
  # x, fixed, model.trend, RMU.data, colname.year=NULL, RMU.names=NULL, index=NULL
  # pt <- list(fixed=result$parametersfixed, RMU.data=result$RMU.data, model.trend=result$model.trend, colname.year=result$colname.year, RMU.names=result$RMU.names)

out <- MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, n.chains = n.chains, n.adapt = n.adapt, 
                 thin=thin, trace=trace, likelihood=fun,
                 fixed=result$parametersfixed, RMU.data=result$RMU.data, model.trend=result$model.trend,
                 colname.year=result$colname.year, RMU.names=result$RMU.names)

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
