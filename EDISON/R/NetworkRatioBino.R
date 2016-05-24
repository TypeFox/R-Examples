#' Calculates the ratio of binomial prior probabilites.
#' 
#' This function calculates the ratio of binomial prior probabilities of two
#' networks.
#' 
#' 
#' @param network.info Network information collected by function
#' \code{\link{CollectNetworkInfo}}.  Note that network.info$new.nets has to be
#' set.
#' @param node.sharing Type of coupling of hyperparameters among nodes:
#' \code{'hard'} or \code{'soft'}.
#' @return Returns the ratio of [prior of new network]/[prior of old network].
#' @author For information about the binomial information sharing prior, see:
#' 
#' Husmeier et al. (2010), "Inter-time segment information sharing for
#' non-homogeneous dynamic Bayesian networks", NIPS.
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @seealso \code{\link{NetworkProbBino}}, \code{\link{CalculatePriorRatio}}
#' @export NetworkRatioBino
NetworkRatioBino <-
function(network.info, node.sharing) {
  # Calculate the ratio of probabilities when applying one edge change to a 
  # network segment
  #
  # Args:
  #   network_info: The network structures and associated information.
  #             network.info$nets         - Structure of all segments 
  #             network.info$prior.params - Hyperparameters alpha, alpha bar, 
  #                                         gamma, gamma bar of the binomial prior
  #             network.info$segment      - Segment being changed
  #             network.info$target       - Target node whose edge is being changed
  #             network.info$parent       - Parent being changed
  #  node.sharing: Indicator flag to decide between 'hard' and 'soft' coupling 
  #                over nodes.
  # Returns:
  #   Ratio of the structure priors

  logprior.old = NetworkProbBino(network.info, 
                     node.sharing)
  
  # Set the network to the new (proposed network). 
  network.info$nets = network.info$new.nets
                   
  logprior.new = NetworkProbBino(network.info, node.sharing)
   
  return(exp(logprior.new - logprior.old));
}

