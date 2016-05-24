#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Arborescence Problems                                               #
#-----------------------------------------------------------------------------#

# getZeroArcs -----------------------------------------------------------------
#' Selects zero weight arcs of a graph
#'
#' Given a directed graph, \code{getZeroArcs} returns the list of arcs with 
#' zero weight. Removes other arcs by assign them infinite value.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return The \code{getZeroArcs} function returns a matrix with the list of
#' zero weight arcs of the graph.
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{msArborEdmonds} and \link{getMinimumArborescence}.

getZeroArcs <- function(nodes, arcs) {
  
  # Work with cost matrix
  Cmat <- ArcList2Cmat(nodes, arcs)
  
  for (i in seq(along = nodes)) {
    # Check all nodes and remove arcs with no zero cost
    Cmat[which(Cmat[, i] != min(Cmat[, i], na.rm = T)), i] <- Inf
  }
  
  # Return list of arcs
  arcs <- Cmat2ArcList(nodes, Cmat, directed = TRUE)
  
  return(arcs)
  
}
#-----------------------------------------------------------------------------#