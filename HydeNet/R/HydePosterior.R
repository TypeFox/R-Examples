#' @name HydePosterior
#' @export HydePosterior
#' 
#' @title Posterior Distributions of a Decision Network
#' @description The posterior distributions of the decision network can be 
#'   evaluated to determine the probabilistic outcomes based on the decision
#'   inputs in the model as well as subject specific factors.
#'   
#' @param cHN A \code{compiledHydeNetwork} object as returned by 
#'   \code{compileJagsNetwork}.
#' @param variable.names a character vector giving the names of variables to be monitored.
#' @param n.iter number of iterations to monitor.
#' @param thin thinning interval for monitors.
#' @param ... options arguments that are passed to the update method for 
#'   jags model objects.
#' @param monitor_observed If TRUE, the observed or fixed variables (those
#'   passed to the \code{data} argument in \code{compileJagsNetwork}) are 
#'   forced into \code{variable.names} if not already provided.  This is 
#'   recommended, especially if you will be binding multiple JAGS runs 
#'   together.
#' @param bind Logical. If \code{TRUE}, posterior distributions will be bound into 
#'   a single data frame.  If \code{FALSE}, the standard output from \code{rjags}
#'   is returned.
#'   
#' @details This is essentially a wrapper around \code{coda.samples} that 
#'   returns in a list the output for each run of \code{coda.samples} over 
#'   the rows of the policy/decision matrix given in the \code{data} argument 
#'   of \code{compileJagsNetwork}.
#'   
#' @return A list of class \code{HydePosterior} with elements \code{codas} 
#'   (the MCMC matrices from \code{coda.samples}), \code{observed} (the values
#'   of the variables that were observed), \code{dag} (the dag object for 
#'   convenience in displaying the network), and \code{factorRef} (giving the
#'   mappings of factor levels to factor variables).  
#'   
#'   The only rationale for giving this object its own class was because it 
#'   produces an enormous amount of material to be printed.  A distinct 
#'   \code{print} method has been written for this object.
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' 
#' @examples
#' data(PE, package="HydeNet")
#' Net <- HydeNetwork(~ wells + 
#'                      pe | wells + 
#'                      d.dimer | pregnant*pe + 
#'                      angio | pe + 
#'                      treat | d.dimer*angio + 
#'                      death | pe*treat,
#'                      data = PE) 
#'   
#'                  
#' compiledNet <- compileJagsModel(Net, n.chains=5)
#' 
#' #* Generate the posterior distribution
#' Posterior <- HydePosterior(compiledNet, 
#'                            variable.names = c("d.dimer", "death"), 
#'                            n.iter = 1000)
#' 
#' #* Posterior Distributions for a Decision Model
#' Net <- setDecisionNodes(Net, angio, treat)
#' decisionNet <- compileDecisionModel(Net, n.chains=5)
#' decisionsPost <- HydePosterior(decisionNet, 
#'                                variable.names = c("d.dimer", "death"),
#'                                n.iter = 1000)
#' 
#' 

HydePosterior <- function(cHN, variable.names, n.iter, thin=1, ...,
                          monitor_observed=TRUE, bind=TRUE){
  if (monitor_observed){
    variable.names <- 
      if (class(cHN$jags) == "jags")
        unique(c(variable.names, names(cHN$observed)))
      else unique(c(variable.names, names(cHN[[1]]$observed)))
  }
    
  
  
  if (class(cHN$jags) == "jags"){
    codas <- coda.samples(cHN$jags, variable.names, n.iter, thin, ...)
    HydePost <- list(codas = codas,
                     observed = cHN$observed,
                     dag = cHN$dag,
                     factorRef = cHN$factorRef)
  }
  else{ 
    codas <- lapply(1:length(cHN),
                    function(j, ...){
                      rjags::coda.samples(cHN[[j]]$jags,
                                    variable.names = variable.names,
                                    n.iter = n.iter,
                                    thin = thin,
                                    ...)
                    },
                    ...)
    observed <- do.call("rbind", lapply(cHN,
                                        function(x) x$observed))
    
    HydePost <- list(codas = codas,
                     observed = observed,
                     dag = cHN[[1]]$dag,
                     factorRef = cHN[[1]]$factorRef)
  }
  
  
  
  class(HydePost) <- "HydePosterior"
  if (bind) return(bindPosterior(HydePost)) else return(HydePost)
  
}
