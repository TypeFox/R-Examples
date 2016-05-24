#' Calculated the global changepoint probabilities.
#' 
#' This function calculates the global probability of a changepoint at each
#' measured timepoint, using the node-specific probabilities.
#' 
#' 
#' @param prob.cps Node-specific changepoint probabilities, a NumNodes by
#' NumTimepoints matrix.
#' @return A matrix of length 1 by NumTimepoints, containing the global
#' changepoint probabilities.
#' @author Frank Dondelinger
#' @seealso \code{\link{calculateCPProbabilities}}
#' @export calculateCPPGlobal
calculateCPPGlobal <-
function(prob.cps) {
  global.cps = matrix(0, 1, dim(prob.cps)[2])
  
  for(node in 1:dim(prob.cps)[1]) {
    global.cps = global.cps + prob.cps[node,] - global.cps * prob.cps[node,]
  }
  
  return(global.cps)
}

