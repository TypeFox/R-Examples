#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Arborescence Problems                                               #
#-----------------------------------------------------------------------------#

# getCheapArcs ----------------------------------------------------------------
#' Substracts the minimum weight of the arcs pointing to each node
#'
#' The \code{getCheapArcs} function substracts to each arc of a given graph the
#' value of the minimum weight of the arcs pointing to the same node.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{getCheapArcs} returns a matrix with a new list of arcs.
#' 
#' @seealso This function is an auxiliar function used in  
#' \link{msArborEdmonds} and \link{getMinimumArborescence}.

getCheapArcs <- function(nodes, arcs) {
  
  # Work with cost matrix
  Cmat <- ArcList2Cmat(nodes, arcs)
  
  for (i in seq(along = nodes)) {
    # Check all nodes
    if (!is.na(all(Cmat[, i] == Inf))) {
      # Only consider arcs with no infinite cost
      min.cost <- Cmat[which(Cmat[, i] == min(Cmat[, i], na.rm = T)), i]
      Cmat[, i] <- Cmat[, i] - min.cost[1]  # substract minimum cost to arcs
    }
  }
  
  # Return list of arcs
  arcs <- Cmat2ArcList(nodes, Cmat, directed = TRUE)
  
  return(arcs)
  
}
#-----------------------------------------------------------------------------#