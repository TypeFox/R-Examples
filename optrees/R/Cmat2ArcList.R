#-----------------------------------------------------------------------------#
# Optrees Package                                                             #
# Auxiliar functions                                                          #
#-----------------------------------------------------------------------------#

# Cmat2ArcList ----------------------------------------------------------------
#' Builds the list of arcs of a graph from its cost matrix
#'
#' The \code{Cmat2ArcList} function builds the list of arcs of a graph from a
#' cost matrix that contains the weights of all the arcs.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param Cmat \eqn{n\times n} matrix that contains the weights or costs of the
#' arcs. Row \eqn{i} and column \eqn{j} represents the endpoints of an arc, and
#' the value of the index \eqn{ij} is its weight or cost. If this value is 
#' \code{NA} or \code{Inf} means that there is no arc \eqn{ij}.
#' @param directed logical value indicating whether the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#'
#' @return \code{Cmat2ArcList} returns a matrix with the list of arcs of the
#' graph. Each row represents one arc. The first two columns contain the two
#' endpoints of each arc and the third column contains their weights.

Cmat2ArcList <- function(nodes, Cmat, directed = TRUE) {
  
  # Initialize
  arcs <- c()  # list of arcs
  
  if (directed) {
    # If the arcs of the graph are directed
    for (i in seq(along = nodes)) {  # check every row
      for (j in seq(along = nodes)) {  # and every column
        if (!is.na(Cmat[i, j]) && Cmat[i, j] != Inf) {
          # Save the arc if exists (no NA nor Inf)
          arcs <- rbind(arcs, c(i, j, Cmat[i, j]))
        }
      }
    }
    
    # Column names if directed graph
    colnames(arcs) <- c("head", "tail", "weight")
    
  } else {
    # If the arcs of the graph are not directed
    for (i in seq(along = nodes)) {  # check every row
      for (j in i:length(nodes)) {  # and every column starting from node i
        if (!is.na(Cmat[i, j]) && Cmat[i, j] != Inf) {
          # Save the arc if exists (no NA nor Inf)
          arcs <- rbind(arcs, c(i, j, Cmat[i, j]))
        }
      }
    }
    
    # Column names if undirected graph
    colnames(arcs) <- c("ept1", "ept2", "weight")
  }
  
  # Return matrix with list of arcs
  return(arcs)
  
}
#-----------------------------------------------------------------------------#