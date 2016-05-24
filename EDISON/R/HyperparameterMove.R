#' Make a hyperparameter move.
#' 
#' This function makes a hyperparameter move for the information sharing prior
#' selected (or no move if no information sharing prior is selected).
#' 
#' 
#' @param method The information sharing method used: \code{'poisson'} for the
#' Poisson prior (no information sharing), \code{'exp_soft'} and
#' \code{'exp_hard'} for the exponential information sharing prior with soft or
#' hard information sharing among nodes, respectively, \code{'bino_soft'} and
#' \code{'bino_hard'} for the binomial information sharing prior with soft or
#' hard information sharing among nodes, respectively.
#' @param network.info Network information collected using
#' \code{\link{CollectNetworkInfo}}.
#' @param GLOBvar Global variables used during the MCMC.
#' @param hyper.proposals Proposal width for hyperparameters.
#' @return List summing up the result of the hypermove. Contains at least:
#' \item{move.made}{Whether a hyperparameter move has been made.}
#' \item{network.info}{The network information, possibly updated if the
#' hyperparameter move was made and accepted.} May contain further elements
#' depending on the type of information sharing prior used. See the
#' prior-specific functions \code{\link{ExpHyperMove}} and
#' \code{\link{BinoHyperMove}} for details.
#' @author Frank Dondelinger
#' @references For information about the information sharing priors, see:
#' 
#' Husmeier et al. (2010), "Inter-time segment information sharing for
#' non-homogeneous dynamic Bayesian networks", NIPS.
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export HyperparameterMove
HyperparameterMove <-
function(method, network.info, GLOBvar, hyper.proposals) {
  # Makes a hyperparameter move (with a certain probability)
  #
  # Args:
  #   method: The information sharing method used. Currently supports "poisson"
  #           (no information sharing) and "exp_soft" (sequential exponential 
  #           information sharing with soft coupling of nodes)
  #   network_info: The network structures and associated information.
  #             network.info$nets        - Structure of all segments 
  #             network.info$target.nets - Structure of all segments in per-node
  #                                        form (easier to use in some
  #                                        situations)
  #             network.info$betas       - Beta parameters for all segments
  #             network.info$segment     - Segment being changed
  #             network.info$target      - Target node whose edge is being 
  #                                        changed
  #             network.info$parent      - Parent being changed
  #
  # Returns:
  #   Structure containing the updated hyperparameter
  
  if(length(network.info$nets) == 1) {
    # Only one segment present, no hyper moves needed
    return(list(move.made=0, network.info)) 
  }

  if(method == 'poisson' || GLOBvar$hyper.fixed) {
    hyper.move = list(move.made=0, network.info)
  } else if(method == 'exp_soft') {
    hyper.move = ExpHyperMove(network.info, 'soft', GLOBvar, hyper.proposals)
  } else if(method == 'exp_hard') {
    hyper.move = ExpHyperMove(network.info, 'hard', GLOBvar, hyper.proposals)
  } else if(method == 'bino_soft') {
    hyper.move = BinoHyperMove(network.info, 'soft', GLOBvar)
  } else if(method == 'bino_hard') {
    hyper.move = BinoHyperMove(network.info, 'hard', GLOBvar)
  } 
  
  return(hyper.move)
}

