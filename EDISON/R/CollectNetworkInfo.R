#' Collects all the network information in one list.
#' 
#' This function collects information about the current network segments and
#' hyperparameters for the information sharing priors.
#' 
#' 
#' @param Sall Structure of all segments. A list of length \code{q}, where each
#' element is a \code{K_i} by \code{q} matrix containing the parents for the
#' current node in each of the \code{K_i} segments.
#' @param Eall Positions of segment boundaries. A list of length \code{q},
#' where each element is a vector containing the segment boundaries for the
#' current parent node.
#' @param prior.params The hyperparameters of the information sharing prior (if
#' applicable).
#' @param posPhase The segment being changed.
#' @param target The target parent node whose edge is being changed.
#' @param q The total number of nodes in the network.
#' @param self.loops Whether self-loops are allowed in the network.
#' @param k The level-2 hyperparameter for the exponential prior.
#' @return The function returns a list with the following elements:
#' \item{nets}{The structure of all segments, a list of length \code{K} where K
#' is the total number of segments over all nodes.} \item{segment}{Identical to
#' \code{posPhase}.} \item{target.nets}{Identical to \code{Sall}.}
#' \item{prior.params}{Identical to \code{prior.params}.}
#' \item{self.loops}{Identical to \code{self.loops}.} \item{k}{Identical to
#' \code{k}.} \item{new.nets}{Dummy variable for holding the proposed network
#' in a network structure move. Originally identical to variable \code{nets}.}
#' @author Frank Dondelinger
#' @export CollectNetworkInfo
CollectNetworkInfo <-
function(Sall, Eall, prior.params, posPhase, target, 
                               q, self.loops, k) {
  # Collects information about the current network segments and hyperparameters
  # for the information sharing priors
  #
  # Args:
  #   Sall:           Structure of all segments
  #   Eall:           Position of segment boundaries
  #   prior.params:   Parameters for the information sharing prior
  #   posPhase:       Segment being changed
  #   target:         Target node whose edge is being changed
  #
  # Returns:
  #   network_info: The network structures and associated information.
  #             network.info$nets         - Structure of all segments
  #             network.info$new.nets     - Structure of all segments after
  #                                         applying the current move 
  #             network.info$target.nets  - Structure of all segments in per-node
  #                                         form (easier to use in some
  #                                         situations)
  #             network.info$prior.params - Parameters for the information
  #                                         sharing prior                                        
  #             network.info$segment      - Segment being changed
  #             network.info$target       - Target node whose edge is being 
  #                                         changed
  #             network.info$self.loops   - Indicator variable showing whether
  #                                         self loops are allowed.
  
  if(any(Sall[[1]][Sall[[1]] != 0] != 1)) 
     stop('Illegal Network Structure in CollectNetworkInfo.') 
  
  network.info = list()
  
  # Convert from list of target nodes to list of segments (easier to 
  # calculate segment similarity)
  converted = convert_nets(Sall, Eall)
  
  network.info$nets = converted$B_nets; 
  network.info$segment = posPhase; network.info$target = target
  network.info$target.nets = Sall
  network.info$prior.params = prior.params
  network.info$self.loops = self.loops
  network.info$k = k
  network.info$cps = Eall
  network.info$global.mapping = converted$mapping
  network.info$global.cps = converted$segs
  network.info$q = q

  for(i in 1:length(network.info$nets)) {
    network.info$nets[[i]] = network.info$nets[[i]][1:q,]
  }
  

  for(j in 1:length(network.info$target.nets)) {
    network.info$target.nets[[j]] = network.info$target.nets[[j]][,1:q,drop=FALSE]
  }

  network.info$new.nets = network.info$nets

  return(network.info)
}

