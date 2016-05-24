#-----------------------------------------------------------------------------#
# Optrees Package                                                             #
# Auxiliar functions                                                          #
#-----------------------------------------------------------------------------#

# removeMultiArcs -------------------------------------------------------------
#' Remove multi-arcs with no minimum cost
#' 
#' The \code{removeMultiArcs} function go through the arcs list of a given 
#' graph and check if there are more than one arc between two nodes. If exist
#' more than one, the function keeps one with minimum cost and remove the 
#' others.
#' 
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param directed logical value indicating whether the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#' 
#' @return \code{removeMultiArcs} returns a new list of arcs without any of
#' the multi-arcs founded.

removeMultiArcs <- function(arcs, directed = TRUE) {
  
  # Intitalize
  new.arcs <- matrix(ncol = 3)[-1, ]  # new list of arcs
  i <- 1  # start with first arc
  
  if (directed) {
    # If the arcs are directed compare only in one direction
    while (nrow(arcs) > 0) {
      # Iterate until all arcs have been checked
      arc <- arcs[i, ]  # select one and check
      if (any(arc[1] == arcs[-i, 1] & arc[2] == arcs[-i, 2])) {
        # If there is more than one arc between two nodes select all of them
        coincidence <- which(arc[1] == arcs[, 1] & arc[2] == arcs[, 2])
        coincidence.arcs <- arcs[coincidence, ]
        
        # Keep the arc with minimum cost
        k <- which(coincidence.arcs[, 3] == min(coincidence.arcs[, 3]))[1]
        min.arc <- coincidence.arcs[k, ]
        new.arcs <- rbind(new.arcs, min.arc)
        # Remove the others
        arcs <- matrix(arcs[-coincidence, ], ncol = 3)
        
      } else {
        # If there is only one arc between two nodes save it directly
        new.arcs <- rbind(new.arcs, arc)
        arcs <- matrix(arcs[-i, ], ncol = 3)  # and remove it from initial list
      }
    }
    
    # Column names if directed
    colnames(new.arcs) <- c("head", "tail", "weight")
    
  } else {
    # If the arcs are undirected compare in two directions
    while (nrow(arcs) > 0) {
      # Iterate until all arcs have been checked
      arc <- arcs[i, ]  # select one and check
      if (any(arc[1] == arcs[-i, 1] & arc[2] == arcs[-i, 2]) | 
            any(arc[1] == arcs[-i, 2] & arc[2] == arcs[-i, 1])) {
        # If there is more than one arc between two nodes select all of them
        coincidence <- which(arc[1] == arcs[, 1] & arc[2] == arcs[, 2] |
                               arc[1] == arcs[, 2] & arc[2] == arcs[, 1])
        coincidence.arcs <- arcs[coincidence, ]
        
        # Keep the arc with minimum cost
        k <- which(coincidence.arcs[, 3] == min(coincidence.arcs[, 3]))[1]
        min.arc <- coincidence.arcs[k, ]
        new.arcs <- rbind(new.arcs, min.arc)
        # Remove the others
        arcs <- matrix(arcs[-coincidence, ], ncol = 3)
        
      } else {
        # If there is only one arc between two nodes save it directly
        new.arcs <- rbind(new.arcs, arc)
        arcs <- matrix(arcs[-i, ], ncol = 3)  # and remove it from initial list
      }
    }
    
    # Column names if undirected
    colnames(new.arcs) <- c("ept1", "ept2", "weight")
    
  }
  
  # Remove row names
  rownames(new.arcs) <- NULL
  
  # Return new list of arcs without multi-arcs
  return(new.arcs)
  
}
#-----------------------------------------------------------------------------#