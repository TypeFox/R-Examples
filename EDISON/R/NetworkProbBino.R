#' Calculates the prior probability of the network segments under the binomial
#' prior.
#' 
#' This function calculates the (log) probability of the network segments using
#' the binomial information sharing prior.
#' 
#' 
#' @param network.info Network information collected by function
#' \code{\link{CollectNetworkInfo}}.
#' @param node.sharing Coupling of hyperparameters among nodes: \code{'hard'}
#' or \code{'soft'}.
#' @return Returns the log prior probability of the network segments under the
#' binomial prior.
#' @author Frank Dondelinger
#' @seealso \code{\link{NetworkRatioBino}}, \code{\link{CalculatePriorRatio}}
#' @references For information about the binomial information sharing prior,
#' see:
#' 
#' Husmeier et al. (2010), "Inter-time segment information sharing for
#' non-homogeneous dynamic Bayesian networks", NIPS.
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export NetworkProbBino
NetworkProbBino <-
function(network.info, node.sharing='soft') {
  # Calculate the (log) probability \prod_i P(M_i|M_{i-1}), where M_i is
  # network segment being changed.
  #
  # Args:
  #   network_info: The network structures and associated information.
  #             network.info$nets         - Structure of all segments 
  #             network.info$prior.params - Hyperparameters alpha, alpha bar, 
  #                                         gamma, gamma bar of the binomial 
  #                                         prior
  #             network.info$segment      - Segment being changed
  #             network.info$target       - Target node whose edge is being changed
  #             network.info$parent       - Parent being changed
  #  node.sharing: Indicator flag to decide between 'hard' and 'soft' coupling 
  #                over nodes.
  # Returns:
  #   P(M_i|M_{i-1})*P(M_{i+1}|M_i)
  
  if(is.null(network.info$prior.params))
    stop(paste("Parameters for prior not set in",
               "NetworkProbBino."))
  
  param.alpha     = network.info$prior.params[1]
  param.alpha.bar = network.info$prior.params[2]
  param.gamma     = network.info$prior.params[3]
  param.gamma.bar = network.info$prior.params[4]
  
  if(node.sharing == 'soft') {
    suff.statistics = CalculateChanges(network.info, node.sharing)
  } else if(node.sharing == 'hard') {
    suff.statistics = CalculateChanges(network.info, node.sharing)
  }
  
  N_1_1 = suff.statistics[1]; N_0_1 = suff.statistics[2]
  N_1_0 = suff.statistics[3]; N_0_0 = suff.statistics[4]
 
  logprior.1 = lgamma(param.alpha + param.alpha.bar) - 
              (lgamma(param.alpha) + lgamma(param.alpha.bar))
  logprior.2 = lgamma(N_1_1 + param.alpha) + lgamma(N_0_1 + param.alpha.bar) -
               lgamma(N_1_1 + param.alpha + N_0_1 + param.alpha.bar)
  logprior.3 = lgamma(param.gamma + param.gamma.bar) - 
              (lgamma(param.gamma) + lgamma(param.gamma.bar))
  logprior.4 = lgamma(N_0_0 + param.gamma) + lgamma(N_1_0 + param.gamma.bar) -
               lgamma(N_0_0 + param.gamma + N_1_0 + param.gamma.bar)            
    
  return(logprior.1 + logprior.2 + logprior.3 + logprior.4)
}

