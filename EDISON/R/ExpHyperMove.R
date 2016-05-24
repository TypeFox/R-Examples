#' Makes an exponential hyperparameter move.
#' 
#' This function tries to make a level-1 or level-2 hyperparameter move for the
#' exponential prior
#' 
#' 
#' @param network.info The network information collected by
#' \code{\link{CollectNetworkInfo}}.
#' @param node.sharing The type of information sharing among nodes:
#' \code{'soft'} or \code{'hard'}.
#' @param GLOBvar Collection of global variables of the MCMC.
#' @param hyper.proposals Proposal width of the hyperparameter move.
#' @return Returns a list with elements: \item{move.made}{1 if a level-1
#' hyperparameter move has been made, 0 otherwise.} \item{network.info}{Network
#' information with updated hyperparameters if the move was accepted.}
#' \item{accept}{Whether a level-1 hyperparameter move has been accepted or
#' not.} \item{move.made.k}{1 if a level-2 hyperparameter move has been made, 0
#' otherwise.} \item{accept.k}{Whether a level-2 hyperparameter move has been
#' accepted or not.} \item{move}{Type of move: 2 for a level-1 hyperparameter
#' move, 3 for a level-2 hyperparameter move.}
#' @author Frank Dondelinger
#' @seealso \code{\link{ExpHyperRatioTarget}}
#' @references For information about the exponential information sharing prior,
#' see:
#' 
#' Husmeier et al. (2010), "Inter-time segment information sharing for
#' non-homogeneous dynamic Bayesian networks", NIPS.
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export ExpHyperMove
ExpHyperMove <-
function(network.info, node.sharing, GLOBvar, hyper.proposals) {
  # Makes a hyperparameter move (with a certain probability)
  #
  # Args:
  #   network_info: The network structures and associated information.
  #             network.info$nets         - Structure of all segments 
  #             network.info$target.nets  - Structure of all segments in per-node
  #                                         form (easier to use in some
  #                                         situations)
  #             network.info$prior.params - Beta parameters for all segments
  #             network.info$segment      - Segment being changed
  #             network.info$target       - Target node whose edge is being 
  #                                        changed
  #             network.info$parent       - Parent being changed
  #
  # Returns:
  #   Structure containing the updated hyperparameter
  
  move.made = 0
  accept = 0
  move.made.k = 0
  accept.k = 0
  
  # Random value for deciding which hyperparameter move to make    
  u = runif(1,0,1)
  
  # Level 1 hyperparameters
  beta.old = network.info$prior.params[network.info$target]
  beta.new = beta.old
  
  # Level 2 hyperparameters
  phi = 0.1; 

  k.old = network.info$k
  k.new = k.old
  move = 5
 
  uniform.range = 10
  proposal.range = hyper.proposals

  # Level-2 hyperparameter move
  if(u > (1-GLOBvar$pp.l2) && node.sharing == 'soft') {
    move.made = 1
    move.made.k = 1 
    move = 6
    
    # Uniform proposal
    k.proposed = proposeContinuous(k.old, proposal.range*10, 
                              uniform.range*10)
    
    betas = network.info$prior.params
    
    log.r = 0 
   
    for(i in 1:length(betas)) {
       log.r = log.r + dgamma(betas[i], k.proposed, scale=phi, log=TRUE) -
                       dgamma(betas[i], k.old, scale=phi, log=TRUE) 
    } 
 
    r = exp(log.r)

    if(runif(1,0,1) < r) {
      k.new = k.proposed  
      accept.k = 1
    }
    
  # Level-2 hyperparameter move for hard coupling
  } else if(u > (1-GLOBvar$pp.l2) && node.sharing == 'hard') {
    move = 6
    move.made = 1
    move.made.k = 1 

    k.proposed = proposeContinuous(k.old, proposal.range*10, uniform.range*10)
    
    r = exp(dgamma(beta.old, k.proposed, scale=phi, log=TRUE) -
                       dgamma(beta.old, k.old, scale=phi, log=TRUE))
 
    if(runif(1,0,1) < r) {
      k.new = k.proposed  
      accept.k = 1
    }
    
  # Level-1 hyperparameter move
  } else if(u > (1-GLOBvar$pp.l1)) {
    move.made = 1
    
    # Uniform proposals
    beta.proposed = proposeContinuous(beta.old, proposal.range, 
                                 uniform.range)
    
    prior.ratio = dgamma(beta.proposed, k.old, scale=phi) / 
        dgamma(beta.old, k.old, scale=phi)
    
    if(node.sharing == 'soft') {
      
      # Get subnetwork corresponding to current target node
      target.net = network.info$target.nets[[network.info$target]]
        
      ratio = ExpHyperRatioTarget(beta.proposed, 
                                   beta.old, target.net, 
                                   network.info$self.loops)

      r = ratio * prior.ratio

    } else if(node.sharing == 'hard') {
      
      r = 1
      
      for(node in 1:length(network.info$target.nets)) {
        
        # Get subnetwork corresponding to current target node
        target.net = network.info$target.nets[[node]]
        ratio.temp = ExpHyperRatioTarget(beta.proposed, beta.old, 
                                             target.net, network.info$self.loops)
     
        r = r * ratio.temp
      }

      r = r * prior.ratio
      
    }
    
    u1 = runif(1,0,1)
 
    if(u1 < r) {
      beta.new = beta.proposed
      accept = 1
    }
    
  }
  
  if(node.sharing == 'soft') {
    network.info$prior.params[network.info$target] = beta.new
    network.info$k = k.new
  } else if(node.sharing == 'hard') {
    network.info$prior.params[] = beta.new
    network.info$k = k.new
  }
  
  return(list(move.made=move.made, network.info=network.info, accept=accept,
              move.made.k=move.made.k, accept.k=accept.k, move=move))
}

