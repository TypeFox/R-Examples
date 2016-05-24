#' Calculates the MH ratio of the binomial prior.
#' 
#' This function calculates the ratio of the binomial information sharing prior
#' with the proposed new hyperparameter values, and the binomial prior with the
#' current hyperparameter values.
#' 
#' 
#' @param params.proposed The new proposed hyperparameter values for the
#' binomial prior.
#' @param changed Gives the index of the parameter that has changed.
#' @param node.sharing Type of information sharing among nodes: \code{'soft'}
#' or \code{'hard'}.
#' @param network.info The network information as collected by
#' \code{\link{CollectNetworkInfo}}.
#' @return This function returns a number greater than zero which represents
#' the ratio of binomial priors.
#' @author Frank Dondelinger
#' @seealso \code{\link{BinoHyperMove}}
#' @references For information about the binomial information sharing prior,
#' see:
#' 
#' Husmeier et al. (2010), "Inter-time segment information sharing for
#' non-homogeneous dynamic Bayesian networks", NIPS.
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export BinoHyperRatio
BinoHyperRatio <-
function(params.proposed, changed, node.sharing,
                                network.info) {
  # Calculates the acceptance ratio of a level-1 hyperparameter move for the 
  # soft binomial prior.
  #
  # Args:
  #  params.proposed: Proposed new hyperparameter value
  #  changed:         Index of the changed parameter
  #  node.sharing:    Indicator flag for soft or hard information sharing among 
  #                   nodes
  #  network_info:    The network structures and associated information.
  #             network.info$nets         - Structure of all segments 
  #             network.info$target.nets  - Structure of all segments in per-node
  #                                         form (easier to use in some
  #                                         situations)
  #             network.info$prior.params - Shared hyperparameters
  #             network.info$segment      - Segment being changed
  #             network.info$target       - Target node whose edge is being 
  #                                        changed
  #             network.info$parent       - Parent being changed
  #
  # Returns:
  #   Acceptance ratio
  
  
  if(node.sharing == 'hard') {
    logprior.old = NetworkProbBino(network.info, node.sharing)
  
    network.info$prior.params[changed] = params.proposed
    
    logprior.new = NetworkProbBino(network.info, node.sharing)
  } else {
    network.info.old = network.info
    network.info.new = network.info
    
    network.info.new$prior.params[changed] = params.proposed
    
    logprior.old = 0
    logprior.new = 0
    
    for(target in 1:dim(network.info$nets[[1]])[1]) {
      network.info.old$target = target
      network.info.new$target = target
      
      logprior.old = logprior.old + NetworkProbBino(network.info.old, 
                                      node.sharing)
      logprior.new = logprior.new + NetworkProbBino(network.info.new, 
                                      node.sharing)
    }
  }
  
  return(exp(logprior.new - logprior.old))
}

