#-----------------------------------------------------------------------------#
# Optrees Package                                                             #
# Auxiliar functions                                                          #
#-----------------------------------------------------------------------------#

# getComponents ---------------------------------------------------------------
#' Connected components of a graph
#'
#' The \code{getComponents} function returns all the connected components of
#' a graph.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{getComponents} returns a list with all the components and the
#' nodes of each one (\code{$components}) and a matrix with all the arcs of the
#' graph and its component (\code{$components.arcs}).
#'  
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,1, 1,6,1, 3,4,1, 4,5,1), ncol = 3, byrow = TRUE)
#' # Components
#' getComponents(nodes, arcs)
#' 
#' @export

getComponents <- function(nodes, arcs) {
  
  # Rebuild graph adding components
  comp.arcs <- matrix(ncol = 4)[-1, ]
  
  # Iterates on each arc starting with first
  i <- 1
  
  # Component nodes
  arc.nodes <- arcs[i, 1:2]
  # List to store component nodes
  comp.nodes <- list(arc.nodes)
  
  # Add first arc and start first component
  comp.arcs <- rbind(comp.arcs, c(arcs[i, ], 1))
  
  # Continue with next arc
  i <- i + 1
  
  # And iterate until check all the arcs
  while (i <= nrow(arcs)) {
    
    # Extract nodes from the arc
    arc.nodes <- arcs[i, 1:2]
    
    # Check nodes and components
    fix.arc <- NULL  # start with arc not fixed to any component
    
    # Iterate for each component
    k <- 1
    nComp <- length(comp.nodes)
    
    while (k <= nComp) {
      
      # Check if there is node in k component
      if (any(arc.nodes %in% comp.nodes[[k]])) {
        
        
        # If arc was already fixed means that it links two components
        if (!is.null(fix.arc)) {
                    
          # Links k component with component saved in fix.arc
          new.nodes <- comp.nodes[[k]]  # retrieves nodes
          # Remove node that already exists in the component
          new.nodes <- new.nodes[-which(new.nodes %in% comp.nodes[[fix.arc]])]
          # Add nodes to corresponding component
          comp.nodes[[fix.arc]] <- c(comp.nodes[[fix.arc]], new.nodes)
          # Removes k component recently linked
          comp.nodes <- comp.nodes[-k]
          
          # If there is one only component makes a new list
          if (class(comp.nodes) == "numeric") {
            comp.nodes <- list(comp.nodes)
          }
          
          # Change number of components in the arcs list
          iArcs <- which(comp.arcs[, 1] %in% comp.nodes[[1]] 
                         & comp.arcs[, 2] %in% comp.nodes[[1]])
          comp.arcs[iArcs, 4] <- fix.arc
          
        } else {
          
          # If there is a node means that the arc belongs to the k component
          comp.arcs <- rbind(comp.arcs, c(arcs[i, ], k))
          # Add new node to the node list of the component
          newNode <- arc.nodes[which(!(arc.nodes %in% comp.nodes[[k]]))]
          comp.nodes[[k]] <- c(comp.nodes[[k]], newNode)
          # Mark arc fixed to component k
          fix.arc <- k
          
        }
        
      }
      
      # Repeat to check if there is another component with the other node
      k <- k + 1
      nComp <- length(comp.nodes)
      
    }
    
    # It there isn't any node in a previous component makes a new one
    if (is.null(fix.arc)) {
      
      # Get number of components and adds one
      comp.nodes[[length(comp.nodes) + 1]] <- arc.nodes
      # Add arc with new component
      comp.arcs <- rbind(comp.arcs, c(arcs[i, ], length(comp.nodes)))
      
    }
    
    # Keep checking until there is no arcs remaining
    i <- i + 1
    
  }
  
  # Set column names in the arcs list
  colnames(comp.arcs) <- c("ept1", "ept2", "weight", "component")
  
  return(list("component" = comp.nodes, "components.arcs" = comp.arcs))
  
}
#-----------------------------------------------------------------------------#