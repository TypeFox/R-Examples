#-----------------------------------------------------------------------------#
# Optrees Package                                                             #
# Auxiliar functions                                                          #
#-----------------------------------------------------------------------------#

# checkGraph ------------------------------------------------------------------
#' Checks if the graph contains at least one tree or one arborescence
#'
#' The \code{checkGraph} function checks if it is possible to find at least one
#' tree (or arborescence, if it is the case) in the graph. It only happens when
#' the graph is connected and it is posible to find a walk from the source to 
#' any other node.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param source.node number pointing the source node of the graph.
#' @param directed logical value indicating whether the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#' 
#' @return \code{checkGraph} returns the value \code{TRUE} if the graph meets
#' the requirements and \code{FALSE} otherwise. If the graph is not acceptable
#' this functions also prints the reason.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,2, 1,3,15, 2,3,1, 2,4,9, 3,4,1), 
#'                byrow = TRUE, ncol = 3)
#' # Check graph
#' checkGraph(nodes, arcs)
#' 
#' @export

checkGraph <- function(nodes, arcs, source.node = 1, directed = TRUE) {
  
  if (is.null(arcs)) {
    # If no arcs graph not valid
    messGraph <- "Empty graph"
    warning(messGraph)
    graphOk <- FALSE
    
  } else if (!all(nodes %in% c(1, arcs[, 1], arcs[, 2])) || 
               !any(arcs[, 1] == source.node)) {
    # If no arcs reaching all nodes, or no arcs leaving source, graph not valid
    messGraph <- "Not connected graph"
    warning(messGraph)
    graphOk <- FALSE
    
  } else {
    # Check walks between nodes
    walks <- c()
    for (i in nodes[-source.node]) {  # i points the final node of the wal
      walk <- searchWalk(nodes, arcs, directed, source.node, i)$walk; walk
      walks <- c(walks, walk)
    }
    
    if (!all(walks)) {
      # If no walk between source node and all of the rest
      messGraph <- "No walk between all the nodes"
      warning(messGraph)
      graphOk <- FALSE
      
    } else {
      # In any other case the graph should be valid
      graphOk <- TRUE
    }
    
  }
  
  return(graphOk)
    
}
#-----------------------------------------------------------------------------#