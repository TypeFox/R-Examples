#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimum Cut Tree Problems                                                   #
#-----------------------------------------------------------------------------#

# searchFlowPath --------------------------------------------------------------
#' Find a maximum flow path
#' 
#' \code{searchFlowPath} go through a given graph and obtains a maximum flow 
#' path between source and sink nodes. The function uses a deep-first search
#' estrategy.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param source.node number pointing to the source node of the graph. It's node
#' \eqn{1} by default.
#' @param sink.node number pointing to the sink node of the graph. It's the
#' last node by default.
#' 
#' @return \code{searchFlowPath} returns a list with two elements:
#' path.nodes vector with nodes of the path.
#' path.arcs matrix with the list of arcs that form the maximum flow path.
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{ghTreeGusfield} and \link{getMinimumCutTree}.

searchFlowPath <- function(nodes, arcs, source.node = 1, 
                           sink.node = nodes[length(nodes)]) {
  
  # Previously add a column with flow equal to the capacity of each arc
  if (ncol(arcs) < 4) {
    arcs <- cbind(arcs, arcs[, 3])
  }
  
  # Initialize
  p.arcs <- matrix(ncol = 4)[-1, ]  # matrix to store arcs of the path
  STOP <- FALSE  # stop condition
  p.nodes <- source.node  # nodes already in the path
  check.nodes <- source.node  # nodes already checked
  
  # Iterate until reach sink node or there isn't path
  while (!(sink.node %in% p.nodes) & !STOP) {
    
    # Select last node in path as actual node
    actual.node <- p.nodes[length(p.nodes)]  # start with source node
    # Not saturated arcs leaving actual node and reaching one not visited
    iArcs <- which(arcs[, 1] == actual.node & arcs[, 4] > 0 
                   & !(arcs[, 2] %in% check.nodes))
    leaving.arcs <- matrix(arcs[iArcs, ], ncol = 4)
    
    if (nrow(leaving.arcs) == 0) {
      # If no leaving arcs go back in path
      p.nodes <- p.nodes[-length(p.nodes)]  # remove last node of the path
      p.arcs <- matrix(p.arcs[-nrow(p.arcs), ], ncol = 4)  # and last arc
      
      if (length(p.nodes) == 0) {
        # If no more nodes there isn't path
        # print("There isn't path")
        STOP <- TRUE
      }
      
    } else {
      # If leaving arcs select one with maximum capacity
      max.arc <- which(leaving.arcs[, 4] == max(leaving.arcs[, 4]))[1]
      # Add arc to path
      p.arcs <- rbind(p.arcs, leaving.arcs[max.arc, ])
      
      # Save reaching node
      p.nodes <- c(p.nodes, leaving.arcs[max.arc, 2])  # in path
      check.nodes <- c(check.nodes, leaving.arcs[max.arc, 2])  # and as checked
    }
    
  }
  
  # Return founded path: nodes and arcs
  output <- list("path.nodes" = p.nodes, "path.arcs" = p.arcs)
  return(output)
  
}
#-----------------------------------------------------------------------------#