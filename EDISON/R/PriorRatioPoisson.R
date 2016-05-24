#' Calculate network prior ratio with Poisson prior.
#' 
#' This function calculates the ratio of the Poisson prior for two networks.
#' 
#' 
#' @param network.info Network information collected using
#' \code{\link{CollectNetworkInfo}}. Note that one needs to set
#' \code{network.info$new.nets}.
#' @param q Number of nodes in the network.
#' @param lambda Vector of lambda hyperparameters for each network.
#' @return Returns the ratio [prior of new network]/[prior of old network].
#' @author Frank Dondelinger
#' @seealso \code{\link{CalculatePriorRatio}}
#' @references For more information on the network structure priors, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export PriorRatioPoisson
PriorRatioPoisson <-
function(network.info, q, lambda) {
  # Calculate the ratio of the Poisson prior for a structure move
  #
  # Args:
  #   network_info: The network structures and associated information.
  #             network.info$nets     - Structure of all segments 
  #             network.info$nets.new - Structure of all segments after making the move
  #             network.info$betas    - Beta parameters for all segments
  #             network.info$segment  - Segment being changed
  #             network.info$target   - Target node whose edge is being changed
  #             network.info$parent   - Parent being changed
  #   q:            Number of nodes in the network
  #   lambda:       Vector of lambda values for each network
  

  ratio = 1

  if(!network.info$self.loops) {
    q = q - 1
  }
  
  # Calculate this prior over local segments only to avoid over-counting
  E = network.info$global.mapping[network.info$target,]
  
  # Calculate for each segment
  for(segment in 1:length(network.info$nets)) {
    
    # Omit repeated segments
    if(segment == 1 || E[segment] != E[segment-1]) {
      s.old = sum(network.info$nets[[segment]][,network.info$target])
      s.new = sum(network.info$new.nets[[segment]][,network.info$target]) 
   
      # If change 
      if(abs(s.old - s.new) != 0) {
        # If added one edge
        if(s.new - s.old == 1) {
          ratio = ratio  * lambda[segment] / (q - s.old)
        # If deleted one edge
        } else if(s.new - s.old == -1) {
          ratio = ratio * (q - s.new) / lambda[segment]
        } else {
          p.new = factorial(q - s.new) * lambda[segment] ^ s.new
          p.old = factorial(q - s.old) * lambda[segment] ^ s.old
          ratio = ratio * (p.new / p.old)
        }
      }  
    }
  }
  return(ratio)
}

