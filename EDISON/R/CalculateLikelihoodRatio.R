#' Calculates the ratio of two likelihoods in a structure move.
#' 
#' This function calculates the ratio of the liklihoods in a network structure
#' move. The returned value is the ratio for the modification of one edge in
#' one segment.
#' 
#' 
#' @param gamma0 Hyperparameter.
#' @param y Target data.
#' @param Pxlm Projection matrix with modified edge.
#' @param Pxl Original projection matrix.
#' @param v0 Hyperparameter.
#' @param delta2 Delta squared parameter (signal-to-noise).
#' @param dir Direction of the change: 1 = Added an edge. 2 = Removed an edge.
#' 0 = No change.
#' @return Returns the likelihood ratio.
#' @author Frank Dondelinger
#' @seealso \code{\link{CalculatePriorRatio}}
#' @references For more information about the hyperparameters and the
#' functional form of the likelihood, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export CalculateLikelihoodRatio
CalculateLikelihoodRatio <-
function(gamma0, y, Pxlm, Pxl, v0, delta2, dir) {
  # Calculate the ratio of the likelihood for a structure move. The returned
  # value is the ratio for the modification of one edge in one segment.
  #
  # Args:
  #   gamma0: Hyperparameter
  #   y: Target data
  #   Pxlm: Projection matrix with modified edge
  #   Pxl: Projection matrix
  #   v0: Hyperparameter
  #   delta2: Current value of delta2 parameter.
  #   dir: Direction of the change (1: Add an edge, -1: Remove an edge, 0: Do nothing)
  #
  # Returns:
  #   Ratio of the likelihoods (for one segment).
  
  r.seg = ((gamma0 + t(y) %*% Pxlm %*% y)/(gamma0 + t(y) %*% Pxl %*% y)) ^
            (-(length(y) + v0)/2)
  
  if(dir != 0) {
    r.seg = r.seg * sqrt(1 + delta2)^(-dir)
  } else {
    r.seg = 1 
  }
          
  return(r.seg)

}

