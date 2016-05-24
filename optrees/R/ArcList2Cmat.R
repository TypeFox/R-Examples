#-----------------------------------------------------------------------------#
# Optrees Package                                                             #
# Auxiliar functions                                                          #
#-----------------------------------------------------------------------------#

# ArcList2Cmat ----------------------------------------------------------------
#' Builds the cost matrix of a graph from its list of arcs
#' 
#' The \code{ArcList2Cmat} function constructs the cost matrix of a graph from
#' a list that contains the arcs and its associated weights.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param directed logical value indicating whether the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#'
#' @return \code{ArcList2Cmat} returns a \eqn{n \times n} matrix that contains
#' the weights of the arcs. It means that the element \eqn{(i,j)} of the matrix
#' returns the weight of the arc \eqn{(i,j)}. If the value of an arc 
#' \eqn{(i,j)} is \code{NA} or \code{Inf}, then it means this arc does not
#' exist in the graph.

ArcList2Cmat <- function(nodes, arcs, directed = TRUE) {
  
  # Initialize
  n <- length(nodes)  # number of nodes
  Cmat <- matrix(nrow = n, ncol = n)  # building cost matrix
  Cmat[, ] <- Inf  # start with all Inf
  
  for (i in 1:nrow(Cmat)) {
    Cmat[i, i] <- NA  # NA loops
  }
  
  for (i in 1:nrow(arcs)) {
    Cmat[arcs[i, 1], arcs[i, 2]] <- arcs[i, 3]  # get cost from list of arcs
  }
  
  if (!directed) {
    # If the arcs of the graph are directed, duplicate arcs
    for (i in 1:nrow(arcs)) {
      Cmat[arcs[i, 2], arcs[i, 1]] <- arcs[i, 3]  # get cost from list of arcs
    }
  }
  
  # Return the cost matrix of the graph
  return(Cmat)
  
}
#-----------------------------------------------------------------------------#