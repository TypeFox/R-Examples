#-----------------------------------------------------------------------------#
# Optrees Package                                                             #
# Auxiliar functions                                                          #
#-----------------------------------------------------------------------------#

# removeLoops -----------------------------------------------------------------
#' Remove loops of a graph
#' 
#' This function reviews the arc list of a given graph and check if exists 
#' loops in it. A loop is an arc that connect a node with itself. If 
#' \code{removeLoops} find a loop remove it from the list of arcs.
#' 
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' 
#' @return \code{removeLoops} returns a new list of arcs without any of the 
#' loops founded.

removeLoops <- function(arcs) {
  
  # Initialize
  i <- 1
  
  while (i <= nrow(arcs)) {
    # Iterate in each row (arcs)
    if (arcs[i, 1] == arcs[i, 2]) {
      # If an arc has same two endpoints remove it
      arcs <- matrix(arcs[-i, ], ncol = 3)
    } else {
      i <- i + 1
    }
    
  }
  
  # Return new list of arcs without loops
  return(arcs)
  
}
#-----------------------------------------------------------------------------#