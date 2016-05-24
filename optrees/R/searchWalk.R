#-----------------------------------------------------------------------------#
# Optrees Package                                                             #
# Auxiliar functions                                                          #
#-----------------------------------------------------------------------------#

# searchWalk ------------------------------------------------------------------
#' Finds an open walk in a graph
#'
#' This function walks a given graph, directed or not, searching for a walk
#' from a starting node to a final node. The \code{searchWalk} function uses a
#' deep-first search strategy to returns the first open walk found, regardless
#' it has formed cycles or repeated nodes.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param directed logical value indicating whether the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#' @param start.node number with one node from which a walk start.
#' @param end.node number with final node of the walk.
#' @param method character string specifying which method use to select the
#' arcs that will form the open walk: \code{"min"} if the function chooses the
#' minimum weight arcs, \code{"max"} if chooses the maximum weight arcs, or 
#' \code{NULL} if chooses the arcs by their order in the list of arcs.
#'
#' @return If \code{searchWalk} found an open walk in the graph returns 
#' \code{TRUE}, a vector with the nodes of the walk and a matrix with the list
#' of arcs of it.

searchWalk <- function(nodes, arcs, directed = TRUE, start.node = nodes[1],
                       end.node = nodes[length(nodes)], method = NULL) {
  
  if (start.node == end.node) {
    stop("Start and end node are the same")
  }
  
  # Previous
  if (!directed) {  # duplicate arcs
    arcs <- rbind(arcs, matrix(c(arcs[, 2], arcs[, 1], arcs[, 3]), ncol = 3))
  }
  # Order the arcs by its endpoints
  arcs <- matrix(arcs[order(arcs[, 1], arcs[, 2]), ], ncol = 3); arcs
  
  # Initialize
  arcs <- cbind(arcs, 0)  # add column to mark checked nodes
  actual.node <- start.node  # start with start.node
  wArcs <- matrix(ncol = 4)[-1, ]  # matrix to store walk arcs
  sWalk <- TRUE  # TRUE to continue searching
  
  while (actual.node != end.node && sWalk) {
    # Keep searching until reach end.node or sWalk == FALSE
    
    # Arcs leaving actual.node
    iArcs <- which(arcs[, 1] == actual.node & arcs[, 4] != 1)  # index
    arcs.leaving <- matrix(arcs[iArcs, ], ncol = 4)  # arcs
    
    if (nrow(arcs.leaving) == 0) {
      # If doesn't exists arcs leaving actual.node back to previous node
      actual.node <- wArcs[nrow(wArcs), 1]  # get previous node
      wArcs <- matrix(wArcs[-nrow(wArcs), ], ncol = 4)  # remove previous arc
      
      if (nrow(wArcs) == 0) {
        # If no arcs in walk stop and get arcs leaving actual.node = source node
        arcs.leaving <- which(arcs[, 1] == actual.node & arcs[, 4] != 1)
        
        if(length(arcs.leaving) == 0) {
          # If no arcs stop
          sWalk <- FALSE  # search ended without walk
        }
      }  
      
    } else {
      # Get arc with choosen method
      if (!is.null(method)) {
        if (method == "max") {  # by maximum weight
          iArc <- which(arcs.leaving[, 3] == max(arcs.leaving[, 3]))[1]
          arc <- arcs.leaving[iArc, ]
        } else if (method == "min") {  # by minimum weight
          iArc <- which(arcs.leaving[, 3] == min(arcs.leaving[, 3]))[1]
          arc <- arcs.leaving[iArc, ]
        } else {
          warning("Unknown method")
        }
      } else {  # by order
        arc <- arcs.leaving[1, ]
      }
      
      # Mark arc as checked and add it to the walk
      arcs[which(arcs[, 1] == arc[1] & arcs[, 2] == arc[2]), 4] <- 1
      wArcs <- rbind(wArcs, arc)
      # Select new actual node
      actual.node <- arc[2]
      
      # Check cycle
      if (actual.node %in% wArcs[, 1]) {
        # Remove last arc
        wArcs <- matrix(wArcs[-nrow(wArcs), ], ncol = 4)
        # Select new actual node
        actual.node <- wArcs[nrow(wArcs), 2]; actual.node
      }
    }
  }
  
  # Get nodes in the founded walk
  wNodes <- unique(c(wArcs[, 1], wArcs[, 2]))
  
  # Prepare list of arcs
  wArcs <- matrix(wArcs[, -4], ncol = 3)  # remove checked column
  rownames(wArcs) <- NULL  # remove row names
  if (sWalk) {
    # If there is a walk change column names
    if (directed) {
      colnames(wArcs) <- c("head", "tail", "weight")  # directed
    } else {
      colnames(wArcs) <- c("ept1", "ept2", "weight")  # not directed
    }
  }
  
  return(list("walk" = sWalk, "walk.nodes" = wNodes, "walk.arcs" = wArcs))
  
}
#-----------------------------------------------------------------------------#