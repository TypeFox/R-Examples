#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimum Cut Tree Problems                                                   #
#-----------------------------------------------------------------------------#

# findstCut -------------------------------------------------------------------
#' Determines the s-t cut of a graph
#' 
#' \code{findstCut} reviews a given graph with a cut between two nodes with the
#' bread-first search strategy and determines the two cut set of the partition.
#' The cut is marked in the arc list with an extra column that indicates the
#' remaining capacity of each arc.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param s number pointing one node of the \code{s} cut in a given graph.
#' It's node \eqn{1} by default.
#' @param t number pointing one node of the \code{t} cut in a given graph.
#' It's the last node by default.
#' 
#' @return \code{findstCut} returns a list with two elements:
#' s.cut vector with nodes of the \code{s} cut.
#' t.cut vector with nodes of the \code{t} cut.
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{ghTreeGusfield} and \link{getMinimumCutTree}.

findstCut <- function(nodes, arcs, s = 1, t = nodes[length(nodes)]) {
  
  if (ncol(arcs) < 4) {
    stop("This function needs an extra column with capacities")
  }
  
  # Remove saturated arcs
  arcs <- matrix(arcs[arcs[, 4] != 0, ], ncol = 4)
  
  # Initialize
  queue.nodes <- s  # start queue with source node
  visit.nodes <- s  # start visited nodes with source node
  
  # Iterate until queue is empty
  while (length(queue.nodes) > 0) {
    
    # Gets first element from queue as actual node
    actual.node <- queue.nodes[1]
    queue.nodes <- queue.nodes[-1]
    
    # Check adjacent nodes
    leaving.arcs <- which(arcs[, 1] == actual.node)  # arcs leaving actual node
    
    if (!(length(leaving.arcs) == 0)) {
      # Make changes only if there are arcs leaving actual node
      adjacent.nodes <- arcs[leaving.arcs, 2]  # adjacent nodes
      
      for (i in 1:length(adjacent.nodes)) {
        # Check every node adjacent to actual node
        new.node <- adjacent.nodes[i]  # select one
        
        if (!(new.node %in% visit.nodes)) {
          # If reaching node not visited select arc
          actualArc <- arcs[which(arcs[, 1] == actual.node 
                                  & arcs[, 2] == new.node), ]
          # Save actual stage
          visit.nodes <- c(visit.nodes, new.node)  # add node to visited
          queue.nodes <- c(queue.nodes, new.node)  # add node to queue
        }
        
      }
      
    }
    
  }
  
  # Recover nodes and build output
  s.cut <- visit.nodes  # reaching nodes from source
  t.cut <- nodes[-which(nodes %in% s.cut)]  # nodes connected to sink
  
  # Return two cuts
  output <- list("s.cut" = s.cut, "t.cut" = t.cut)
  return(output)
  
}
#-----------------------------------------------------------------------------#