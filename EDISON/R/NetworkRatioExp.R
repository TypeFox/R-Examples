#' Calculates the ratio of exponential network prior probabilities.
#' 
#' This function calculates the ratio of exponential network information
#' sharing prior probabilities.
#' 
#' 
#' @param network.info Network information collected using the function
#' \code{\link{CollectNetworkInfo}}. Note that \code{network.info$new.nets} has
#' to be set.
#' @return Returns the ratio [prior of new network]/[prior of old network].
#' @author Frank Dondelinger
#' @seealso \code{\link{NetworkProbExp}}, \code{\link{CalculatePriorRatio}}
#' @references For information about the exponential information sharing prior,
#' see:
#' 
#' Husmeier et al. (2010), "Inter-time segment information sharing for
#' non-homogeneous dynamic Bayesian networks", NIPS.
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export NetworkRatioExp
NetworkRatioExp <-
function(network.info) {
  # Calculate the ratio of probabilities when applying one edge change to a 
  # network segment
  #
  # Args:
  #   network_info: The network structures and associated information.
  #             network.info$nets         - Structure of all segments 
  #             network.info$prior.params - Beta parameters for all segments
  #             network.info$segment      - Segment being changed
  #             network.info$target       - Target node whose edge is being changed
  #             network.info$parent       - Parent being changed
  # Returns:
  #   Ratio of the structure priors

  logprior.old = NetworkProbExp(network.info)
  
  network.info$nets = network.info$new.nets
    
  logprior.new = NetworkProbExp(network.info)
  
  return(exp(logprior.new - logprior.old));
}

