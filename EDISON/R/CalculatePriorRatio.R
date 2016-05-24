#' Calculates the network prior ratio.
#' 
#' This function calculates the ratio of the network structure priors for a
#' structure move.
#' 
#' 
#' @param method Indicates which prior to use: \code{'poisson'} for the
#' standard Poisson prior (no information sharing), \code{'exp_soft'} or
#' \code{'exp_hard'} for the exponential information sharing prior with soft or
#' hard sharing among nodes and \code{'bino_soft'} or \code{'bino_hard'} for
#' the binomial information sharing prior with soft or hard sharing among
#' nodes.
#' @param q Number of nodes in the network.
#' @param lambda Vector of lambda hyperparameter values for each network
#' (needed for the Poisson prior).
#' @param network.info The network information collected using
#' \code{\link{CollectNetworkInfo}}.
#' @return Returns the ratio of the network structure priors for the proposed
#' structure move.
#' @author Frank Dondelinger
#' @seealso \code{\link{CalculateLikelihoodRatio}}
#' @references For more information on the network structure priors, see:
#' 
#' Husmeier et al. (2010), "Inter-time segment information sharing for
#' non-homogeneous dynamic Bayesian networks", NIPS.
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export CalculatePriorRatio
CalculatePriorRatio <-
function(method, q, lambda, network.info) {
  # Calculate the ratio of the network structure priors for a structure move.
  #
  # Args:
  #   method: String showing which prior to use:
  #             "poisson"  - Standard poisson prior (no information sharing)
  #             "exp_soft" - Exponential sequential prior with soft information
  #                          sharingi
  #   q:      Number of nodes in the network
  #   lambda: Vector of lambda values for each network (needed for Poisson prior)
  #   network_info: The network structures and associated information.
  #             network.info$nets     - Structure of all segments 
  #             network.info$betas    - Beta parameters for all segments
  #             network.info$segment  - Segment being changed
  #             network.info$target   - Target node whose edge is being changed
  #             network.info$parent   - Parent being changed
  # Returns:            
  #   Ratio of the structure priors.
  
  if(method == 'poisson') {
    prior.ratio = PriorRatioPoisson(network.info, q, lambda)
  } else if(method == 'exp_soft' || method == 'exp_hard') {
    prior.ratio = NetworkRatioExp(network.info)
  } else if(method == 'bino_soft') {
    prior.ratio = NetworkRatioBino(network.info, 'soft')
  } else if(method == 'bino_hard') {
    prior.ratio = NetworkRatioBino(network.info, 'hard')
  } 
  
  return(prior.ratio)
}

