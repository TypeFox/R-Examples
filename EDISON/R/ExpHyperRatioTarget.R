#' Calculates the ratio of an exponential hyperparameter move.
#' 
#' This function calculates the acceptance ratio of a level-1 hyperparameter
#' move for a given target node.
#' 
#' 
#' @param beta.proposed Proposed new hyperparameter value.
#' @param beta.old Previous value of hyperparameter beta.
#' @param target.net Network segments for the target node associated with this
#' hyperparameter value.
#' @param self.loops \code{'TRUE'} if self-loops are acceptable, \code{'FALSE'}
#' otherwise.
#' @return Returns the ratio of the exponential prior with the previous
#' hyperparameter value and the proposed new hyperparameter value.
#' @author Frank Dondelinger
#' @seealso \code{\link{ExpHyperMove}}
#' @references For information about the exponential information sharing prior,
#' see:
#' 
#' Husmeier et al. (2010), "Inter-time segment information sharing for
#' non-homogeneous dynamic Bayesian networks", NIPS.
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export ExpHyperRatioTarget
ExpHyperRatioTarget <-
function(beta.proposed, beta.old, target.net, self.loops) {
  # Calculates the acceptance ratio of a level-1 hyperparameter move for the 
  # soft exponential prior.
  #
  # Args:
  #  beta.proposed: Proposed new hyperparameter value
  #  beta.old:      Old hyperparameter value
  #  target.net:    Network segments for the target associated with the hyper-
  #                 parameter value
  #
  # Returns:
  #   Acceptance ratio
  
  n.segs = dim(target.net)[1]
  n.parents = dim(target.net)[2]

  if(!self.loops) {
    n.parents = n.parents - 1
  }
   
  # If there is only one segment, there is nothing to calculate
  if(is.null(n.segs)) return(1) 
  
  # Network segments h+1
  target.net.shifted = target.net[2:n.segs,] 
  
  # Calculate differences
  differences = sum(c(abs(target.net.shifted - target.net[1:(n.segs-1),])))
 
  log.ratio.diff = (beta.old - beta.proposed) * differences
 
  # Calculate Z
  log.ratio.z = (n.segs-1)*n.parents*(log(1 + exp(-beta.old)) - 
    log(1 + exp(-beta.proposed)))
     
  return(exp(log.ratio.diff + log.ratio.z))
}

