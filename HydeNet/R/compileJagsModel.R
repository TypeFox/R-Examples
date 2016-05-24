#' @name compileJagsModel
#' @export compileJagsModel
#' @import rjags
#' 
#' @title Compile Jags Model from a Hyde Network
#' @description Generates the JAGS code from the Hyde network and uses it to 
#'   create an object representing a Bayesian graphical model.
#'   
#' @param network An object of class \code{HydeNetwork}
#' @param data A list of data values to be observed in the nodes.  It is
#'   passed to the \code{data} argument of \code{rjags::jags}.  Alternatively,
#'   a data frame representing a policy matrix may be provided to compile
#'   multiple JAGS models.
#' @param ... Additional arguments to be passed to \code{jags.model}
#' 
#' @details \code{compileJagsModel} is a partial wrapper for 
#'   \code{jags.model}. Running \code{compileJagsModel(network)} is 
#'   equivalent to running \code{jags.model(textConnection(writeNetworkModel(network)))}.
#'   
#' @return Returns a \code{compiledHydeNetwork} object.  The \code{jags} element
#'   of this object is suitable to pass to \code{coda.samples}.  Otherwise, 
#'   the primary function of the object is plotting the network with 
#'   observed data shown.   
#'   
#' @author Benjamin Nutter
#' @seealso \code{jags.model} 
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
#' Posterior
#' 
#' #* For a single model (ie, not a decision model), the user may choose to 
#' #* use the \code{rjags} function \code{coda.samples}.
#' #* However, this does not have a succinct print method
#' s <- coda.samples(compiledNet$jags, 
#'                   variable.names = c("d.dimer", "death"), 
#'                   n.iter=1000)
#'                 

compileJagsModel <- function(network, data=NULL, ...){
  
  factorRef <- makeFactorRef(network)
  
  #* convert label to value
  data <- convertLabelToValue(data, factorRef)

  cpt_arrays <- makeCptArrays(network) #* The utilty function is in the 
                                       #* file for compileDecisionModel
 
  jags <- rjags::jags.model(textConnection(writeNetworkModel(network)), 
                    data = if (is.null(data) & length(cpt_arrays) == 0) sys.frame(sys.parent()) 
                             else c(data, cpt_arrays), ...)
  
  #* cHN for compiled Hyde Network
  cHN <- list(jags=jags, observed=data, dag=network$dag, factorRef=factorRef)
  
  class(cHN) <- c("compiledHydeNetwork")
  cHN
}




#****** UTILITY FUNCTIONS
convertLabelToValue <- function(data, factorRef){
  msg <- ""
  for (i in names(data)){
    if (!is.numeric(data[[i]]))
    {
      if (!i %in% names(factorRef)){
        msg <- c(msg,
                 paste0("'", i, "' was not numeric and no matching factor could be found in the data."))
      }
      else if (!all(data[[i]] %in% factorRef[[i]]$label))
      {
        msg <- c(msg,
                 paste0("Values observed in '", i, "' must be one of ",
                        paste0(c(factorRef[[i]]$label, factorRef[[i]]$value), collapse = ", ")))
      }
      else data[[i]] <- factorRef[[i]]$value[which(factorRef[[i]]$label == data[[i]])]
    }
  }

  if (length(msg) > 1) stop(paste(msg, collapse="\n"))
  return(data)
}