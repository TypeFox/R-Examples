#' Function to calculate the number of differences between adjaccent network
#' segments.
#' 
#' This function takes the current network structure, compares each segment to
#' the next one, and calculates the number of changes. If soft information
#' sharing among nodes is active, then this procedure is only done for the
#' current target node.
#' 
#' 
#' @param network.info The network information collected by function
#' \code{\link{CollectNetworkInfo}}.
#' @param node.sharing Specifies the type of information sharing among nodes:
#' \code{'soft'} or \code{'hard'}.
#' @return Returns a vector with 4 elements: the number of coinciding edges,
#' the number of edges in the previous segment that are absent in the next one,
#' the number of edges in the next segment that are absent in the previous one
#' and the number of coinciding non-edges.
#' @author Frank Dondelinger
#' @export CalculateChanges
CalculateChanges <-
function(network.info, node.sharing) {
  # Utility function to calculate the number of differences between adjacent 
  # network segments.
  #
  # Args:
  #   network_info: The network structures and associated information.
  #             network.info$nets         - Structure of all segments 
  #             network.info$prior.params - Parameters alpha, alpha bar, gamma,
  #                                         gamma bar of the binomial prior
  #             network.info$segment      - Segment being changed
  #             network.info$target       - Target node whose edge is being changed
  #             network.info$parent       - Parent being changed
  # Returns:
  #   Matrix with 4 entries corresponding to N_1_1, N_0_1, N_1_0 and N_0_0.
  
  # Initialise sufficient statistics
  N_1_1 = 0; N_0_1 = 0
  N_1_0 = 0; N_0_0 = 0
  
  N_0_0_inactive = 0
  target = network.info$target
  q = length(network.info$target.nets)
   
  if(length(network.info$nets) > 1) {
    for(seg in 2:length(network.info$nets)) {
      # Soft information sharing
      if(node.sharing == 'soft') {
        network.prev = network.info$nets[[seg-1]][,target]
        network.next = network.info$nets[[seg]][,target] 

        N_0_0_inactive = N_0_0_inactive + 1
      # Hard information sharing
      } else if(node.sharing == 'hard') {
        network.prev = network.info$nets[[seg-1]]
        network.next = network.info$nets[[seg]]
        N_0_0_inactive = N_0_0_inactive + q
      }
       
      differences = c(network.prev - 2*network.next)
      
      N_1_1 = N_1_1 + sum(differences == -1)
      N_0_1 = N_0_1 + sum(differences == 1)
      N_1_0 = N_1_0 + sum(differences == -2)
      N_0_0 = N_0_0 + sum(differences == 0)
      
    } 
  }
  
  if(!network.info$self.loops) {
    N_0_0 = N_0_0 - N_0_0_inactive
  }

  result = matrix(0, 4, 1)
  
  result[1] = N_1_1; result[2] = N_0_1
  result[3] = N_1_0; result[4] = N_0_0
  
  return(result)
}

