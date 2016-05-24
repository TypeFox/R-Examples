#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Arborescence Problems                                               #
#-----------------------------------------------------------------------------#

# getMinCostArcs --------------------------------------------------------------
#' Selects the minimum cost of the arcs pointing to each node
#'
#' Given a directed graph, \code{getMinCostArcs} selects the minimum cost arcs
#' entering each node and removes the others.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return The \code{getMinCostArcs} function returns a matrix with the list of
#' the minimum cost arcs pointing to each node of the graph.
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{msArborEdmonds} and \link{getMinimumArborescence}.

getMinCostArcs <- function(nodes, arcs) {
  
  # Work with cost matrix
  Cmat <- ArcList2Cmat(nodes, arcs)
  
  for (i in seq(along = nodes)) {
    # Remove no minimum cost edges assigning infinite cost
    Cmat[which(Cmat[, i] != min(Cmat[, i], na.rm = T)), i] <- Inf
  }
  
  # Return list of arcs
  arcs <- Cmat2ArcList(nodes, Cmat, directed = TRUE)
  
  return(arcs)
  
}
#-----------------------------------------------------------------------------#