#' @name print.HydePosterior
#' @export 
#' @method print HydePosterior
#' 
#' @title Print a Hyde Posterior Distribution Object
#' @details Prints a brief description of a HydePosterior object.
#' 
#' @param x a \code{HydePosterior} object
#' @param ... additional arguments to be passed to print methods.  Currently 
#'   none in use.
#'   
#' @details Prints the number of posterior distributions, chains, and 
#'   iterations, as well as the observed values.
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
#' Net <- setDecisionNodes(Net, treat)  
#'                  
#' compiledNet <- compileJagsModel(Net, n.chains=5)
#' 
#' #* Generate the posterior distribution for the model (but not the 
#' #* decision model)
#' Posterior <- HydePosterior(compiledNet, 
#'                            variable.names = c("d.dimer", "death"), 
#'                            n.iter = 1000)
#' Posterior
#' 
#' #* Generate the posterior for the decision model
#' Decision <- compileDecisionModel(Net, n.chains=5)
#' Posterior_decision <- HydePosterior(Decision, 
#'                                     variable.names = c("d.dimer", "death"), 
#'                                     n.iter = 1000)
#' 

print.HydePosterior <- function(x, ...){
  n_distributions <- if (class(x$codas) == "mcmc.list") 1 else length(x$codas)
  n_chains <- if (class(x$codas) == "mcmc.list") length(x$codas) else length(x$codas[[1]])
  n_iterations <- if (class(x$codas) == "mcmc.list") nrow(x$codas[[1]]) else nrow(x$codas[[1]][[1]])
  
  cat(paste0("Posterior distributions of a Hyde Network\n",
             "number of posterior distributions: ", n_distributions, "\n",
             "number of chains: ", n_chains, "\n",
             "number of iterations: ", n_iterations, "\n"))
  
  if (!is.null(x$observed)){
    cat("\nObserved at the values:\n")
    print(x$observed)
  }
}
