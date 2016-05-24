#' STRN_MHmcmc runs the Metropolis-Hastings algorithm for STRN (Bayesian MCMC)
#' @title Metropolis-Hastings algorithm for Sexualisation Thermal Reaction Norm
#' @author Marc Girondot
#' @return A list with resultMCMC being mcmc.list object, resultLnL being likelihoods and parametersMCMC being the parameters used
#' @param n.iter Number of iterations for each step
#' @param parametersMCMC A set of parameters used as initial point for searching with information on priors
#' @param result An object obtained after a STRN fit
#' @param n.chains Number of replicates
#' @param n.adapt Number of iterations before to store outputs
#' @param thin Number of iterations between each stored output
#' @param trace True or False, shows progress
#' @param batchSize Number of observations to include in each batch fo SE estimation
#' @param dataSTRN A named list data used to estimate likelihoods (see further in description)
#' @param filename If intermediate is not NULL, save intermediate result in this file
#' @param intermediate Period for saving intermediate result, NULL for no save
#' @param previous Previous result to be continued. Can be the filename in which intermediate results are saved.
#' @description Run the Metropolis-Hastings algorithm for Sexualisation Thermal Reaction Norm.\cr
#' The number of iterations is n.iter+n.adapt+1 because the initial likelihood is also displayed.\cr
#' I recommend that thin=1 because the method to estimate SE uses resampling.\cr
#' If initial point is maximum likelihood, n.adapt = 0 is a good solution.\cr
#' To get the SE of the point estimates from \code{result_mcmc <- STRN_MHmcmc(result=try)}, use:\cr
#' \code{result_mcmc$SD}\cr
#' coda package is necessary for this function.\cr
#' The dataSTRN is a named list with the following objects:\cr
#' \itemize{
#'   \item EmbryoGrowthTRN= result of \code{\link{searchR}}
#'   \item tsd= result of \code{\link{tsd}}
#'   \item sexed= vector with number of sexed embryos
#'   \item males= vector with number of males (could be also females=)
#'   \item Temperatures= a text of the temperatures name used as CTE
#' }
#' The Temperatures text for CTE can be:
#' \itemize{
#'   \item \code{TimeWeighted.temperature.mean}
#'   \item \code{TSP.TimeWeighted.temperature.mean}
#'   \item \code{TSP.MassWeighted.temperature.mean}
#'   \item \code{TSP.STRNWeighted.temperature.mean}
#'   \item \code{TSP.MassWeighted.STRNWeighted.temperature.mean}
#'   \item \code{MiddleThird.TimeWeighted.temperature.mean}
#' }
#' They are explained in the \code{\link{info.nests}} function.\cr
#' This function is not still fully described as it has not been published still.\cr
#' The parameters intermediate and filename are used to save intermediate results every 'intermediate' iterations (for example 1000). Results are saved in a file of name filename.\cr
#' The parameter previous is used to indicate the list that has been save using the parameters intermediate and filename. It permits to continue a mcmc search.\cr
#' These options are used to prevent the consequences of computer crash or if the run is very very long and processes at time limited.\cr
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' }
#' @export

STRN_MHmcmc <- function(result=NULL, n.iter=10000, 
parametersMCMC=NULL, n.chains = 1, n.adapt = 0, 
thin=1, trace=FALSE, batchSize=sqrt(n.iter), 
dataSTRN=NULL, 
intermediate=NULL, filename="intermediate.Rdata", previous=NULL) {
  
  if (is.character(previous)) {
    itr <- NULL
    load(previous)
    previous <- itr
    rm(itr)
    print("Continue previous mcmc run")
  } else {
    print(parametersMCMC)
  }

# 29/1/2014; Ajout de result$weight
out <- MHalgoGen(n.iter=n.iter, parameters=parametersMCMC, 
n.chains = n.chains, n.adapt = n.adapt, thin=thin, trace=trace, 
	data=dataSTRN, likelihood=get(".fonctionSTRNMCMC"), 
intermediate=intermediate, filename=filename, previous=previous)

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
