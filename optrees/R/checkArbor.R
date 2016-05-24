#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Arborescence Problems                                               #
#-----------------------------------------------------------------------------#

# checkArbor ------------------------------------------------------------------
#' Checks if there is at least one arborescence in the graph
#'
#' Given a directed graph, \code{checkArbor} searchs for an arborescence from
#' the list of arcs. An arborescence is a directed graph with a source node and
#' such that there is a unique path from the source to any other node.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param source.node source node of the graph. Its default value is \eqn{1}.
#'
#' @return If \code{checkArbor} found an arborescence it returns \code{TRUE}, 
#' otherwise it returns \code{FALSE}. If there is an arborescence the function 
#' also returns the list of arcs of the arborescence.
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{msArborEdmonds} and \link{getMinimumArborescence}.

checkArbor <- function(nodes, arcs, source.node = 1) {
  
  # Previous
  arcs <- cbind(arcs, 0)  # add column to mark arcs checked
  arbor.nodes <- c()  # vector to store nodes of the arboresecence
  walk.nodes <- c()  # vector to store nodes of the actual walk
  arborescence <- FALSE  # start without arborescence
  arbor.arcs <-  matrix(ncol = 4)[-1, ]  # matrix to store arcs of arborescence
  
  # Initialize
  actual.node <- source.node  # always start in source node
  arbor.nodes <- c(arbor.nodes, actual.node)  # add node to arborescence
  walk.nodes <- c(walk.nodes, actual.node)  # add node to walk
  
  # Iterate until find arborescence or remove all nodes from walk
  while (arborescence == FALSE && length(walk.nodes) > 0) {
    
    # Select arcs leaving the current node
    i <- which(arcs[, 1] == actual.node & arcs[, 4] == 0)
    arcsLeaving <- matrix(arcs[i, ], ncol = 4)
    
    if (nrow(arcsLeaving) == 0) {
      # If no arcs leaving the current node go back in walk
      walk.nodes <- walk.nodes[-length(walk.nodes)]  # remove node from walk
      actual.node <- walk.nodes[length(walk.nodes)]  # get previous node
      
    } else {
      # If there are arcs leaving the current node select one
      arc <- matrix(arcsLeaving[1, ], ncol = 4)  # get first arc
      
      # Check the validity of the arc
      if (arc[2] %in% arbor.nodes || arc[2] %in% walk.nodes) {
        # If the arc reaches a node that is already in arborescence or in walk
        # Mark the arc as checked and returns
        arcs[which(arcs[, 1] == arc[1] & arcs[, 2] == arc[2]), 4] <- 1
        # TO-DO: Cycle if the node is in the walk
        
      } else {
        # The arc can be added to the arborescence
        arbor.arcs <- rbind(arbor.arcs, arc)
        # Mark the arc as checked and update nodes
        arcs[which(arcs[, 1] == arc[1] & arcs[, 2] == arc[2]), 4] <- 1
        actual.node <- arc[2]  # reaching node is new actual node
        arbor.nodes <- c(arbor.nodes, actual.node)  # add node to arborescence
        walk.nodes <- c(walk.nodes, actual.node)  # add node to walk
      }
    }
    
    if (all(nodes %in% arbor.nodes)) {
      # If all nodes are in arbor.nodes we have found the arborescence
      arborescence <- TRUE
    }
    
  }
  
  if (arborescence == FALSE) {
    # If no arborescence clear all stored arcs
    arbor.arcs <- matrix(ncol = 4)[-1, ]; arbor.arcs
  }
  
  # Remove last column
  arbor.arcs <- arbor.arcs[, -4]
  
  return(list("arborescence" = arborescence, "arbor.arcs" = arbor.arcs))
  
}
#-----------------------------------------------------------------------------#