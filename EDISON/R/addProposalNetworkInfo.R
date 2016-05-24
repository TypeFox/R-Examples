#' Add the proposed new network to the new.nets list.
#' 
#' Updates the network.info data structure so that it stays consistent.
#' 
#' @param network.info Data structure containing the current network.
#' @param newS Proposed new network for this target, a num.local.segs by 
#' num.parents matrix.
#' @param E The current vector of local segments for this target (only used to 
#' check for consistency with the network.info change points).
#' @return Updated network.info data structure, with new network added to 
#' new.nets.
#' @author Frank Dondelinger
#' @export addProposalNetworkInfo
addProposalNetworkInfo <- function(network.info, newS, E) {
  # Update network info with proposed new network (assumes no changes to the changepoints)
  
  if(any(E != network.info$cps[[network.info$target]])) {
    network.info = NULL
  } else {
    
    global.update = newS[network.info$global.mapping[network.info$target,],,drop=FALSE]
    
    network.info$new.nets = lapply(1:length(network.info$nets),
      function(i) {
        net.temp = network.info$nets[[i]]
        net.temp[,network.info$target] = global.update[i,]
        return(net.temp)
      })
  }
  
  return(network.info)
}