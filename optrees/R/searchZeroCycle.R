#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Arborescence Problems                                               #
#-----------------------------------------------------------------------------#

# searchZeroCycle -------------------------------------------------------------
#' Zero weight cycle in a graph
#'
#' Given a directed graph, \code{searchZeroCycle} search paths in it that forms
#' a zero weight cycle. The function finishes when found one cycle.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{searchZeroCycle} returns a vector with the nodes and a matrix
#' with a list of arcs of the cycle found.
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{msArborEdmonds} and \link{getMinimumArborescence}.

searchZeroCycle <- function(nodes, arcs) {
  
  # Previous
  arcs <- cbind(arcs, 0)  # add column to mark arcs checked
  uncheck.nodes <- nodes  # vector with uncheked nodes
  stack.nodes <- c()  # vector to store nodes of the actual walk
  cycle.nodes <- c()  # vector to store nodes of the cycle
  
  # Initialize
  arc <- arcs[1,]  # start with first arc
  arcs[1, 4] <- 1  # mark arc as checked
  
  # Start with first node of the arc, remove from unchecked and add to walk
  actual.node <- arc[1]
  uncheck.nodes <- uncheck.nodes[-which(uncheck.nodes == actual.node)]
  stack.nodes <- c(stack.nodes, actual.node)
  # Continue with second node of the arc, remove from unchecked and add to walk
  actual.node <- arc[2]
  uncheck.nodes <- uncheck.nodes[-which(uncheck.nodes == actual.node)]
  stack.nodes <- c(stack.nodes, actual.node)
  
  # Iterate until found a cycle and there are unchecked arcs
  while (length(cycle.nodes) == 0 && any(arcs[, 4] != 1)) {
    
    # Check uncheked arcs leaving actual node
    i <- which(arcs[, 1] == actual.node & arcs[, 4] == 0)
    leaving.arcs <- matrix(arcs[i, ], ncol = 4)  # arcs leaving actual node
    # WARNING! The arc may have been checked but reach a valid node
    
    if (nrow(leaving.arcs) == 0) {
      # If no arcs leaving go back in stack
      stack.nodes <- stack.nodes[-length(stack.nodes)]  # remove node
      
      # Review stack
      if (length(stack.nodes) == 0) {
        # If stack is empty find another uncheked arc
        uncheck.arcs <- arcs[which(arcs[, 4] == 0),]
        arc <- matrix(uncheck.arcs[1, ], ncol = 4)  # get first unchecked arc
        arcs[which(arcs[, 1] == arc[1] & arcs[, 2] == arc[2]), 4] <- 1  # mark
        
        # Get first node of the arc, remove from unchecked and add to walk
        actual.node <- arc[1]
        uncheck.nodes <- uncheck.nodes[-which(uncheck.nodes == actual.node)]
        stack.nodes <- c(stack.nodes, actual.node)
        # Get second node of the arc, remove from unchecked and add to walk
        actual.node <- arc[2]
        uncheck.nodes <- uncheck.nodes[-which(uncheck.nodes == actual.node)]
        stack.nodes <- c(stack.nodes, actual.node)
        
      } else {
        # There are nodes in stack so check precedent
        actual.node <- stack.nodes[length(stack.nodes)]
      }
      
    } else {
      # If arcs leaving get one and check cycle
      arc <- matrix(leaving.arcs[1, ], ncol = 4)  # get first arc
      arcs[which(arcs[, 1] == arc[1] & arcs[, 2] == arc[2]), 4] <- 1  # mark
      
      # Check cycle
      if (arc[2] %in% stack.nodes) {
        # If reaching node is in stack there is a cycle
        from <- which(stack.nodes == arc[2])  # stack index where cycle begins
        to <- which(stack.nodes == arc[1])  # stack index where cycle ends
        # Save cycle
        cycle.nodes <- c(cycle.nodes, stack.nodes[from:to], stack.nodes[from])
        actual.node <- arc[1]; actual.node  # keep actual node
        
      } else {  
        # If reaching node isn't in stack no cycle
        actual.node <- arc[2]  # reaching node is the new node to check
        # Remove node from unchecked and add to walk
        uncheck.nodes <- uncheck.nodes[-which(uncheck.nodes == actual.node)]
        stack.nodes <- c(stack.nodes, actual.node)
        
      }
    }
  }
  
  # Final check
  if (length(cycle.nodes) == 0) {
    # If cycle.nodes is empty there is no cycle
    print("No cycle!")
    cycle <- FALSE
    return(cycle)
    
  } else {
    # There is a cycle
    cycle <- TRUE
    cycle.arcs <- matrix(ncol = 4)[-1, ]  # matrix to store arcs of the cycle
    for (i in 1:(length(cycle.nodes)-1)) {
      # Check all nodes in the cycle to add arcs
      i <- which(arcs[, 1] == cycle.nodes[i] & arcs[,2] == cycle.nodes[i+1])
      cycle.arcs <- rbind(cycle.arcs, arcs[i, ])
    }
    
    cycle.arcs <- cycle.arcs[, -4]  # remove last column
    
    return(list("cycle" = cycle, "cycle.nodes" = cycle.nodes,
                "cycle.arcs" = cycle.arcs))
    
  }
}
#-----------------------------------------------------------------------------#