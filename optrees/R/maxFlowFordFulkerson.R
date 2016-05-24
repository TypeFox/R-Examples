#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimum Cut Tree Problems                                                   #
#-----------------------------------------------------------------------------#

# maxFlowFordFulkerson --------------------------------------------------------
#' Maximum flow with the Ford-Fulkerson algorithm
#' 
#' The \code{maxFlowFordFulkerson} function computes the maximum flow in a 
#' given flow network with the Ford-Fulkerson algorithm.
#' 
#' @details The Ford-Fulkerson algorithm was published in 1956 by L. R. Ford, 
#' Jr. and D. R. Fulkerson. This algorithm can compute the maximum flow between
#' source and sink nodes of a flow network.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param directed logical value indicating wheter the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#' @param source.node number pointing to the source node of the graph. It's node
#' \eqn{1} by default.
#' @param sink.node number pointing to the sink node of the graph. It's the
#' last node by default.
#' 
#' @return \code{maxFlowFordFulkerson} returns a list with:
#' s.cut vector with nodes of the \code{s} cut.
#' t.cut vector with nodes of the \code{t} cut.
#' max.flow value with the maximum flow in the flow network.
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{ghTreeGusfield} and \link{getMinimumCutTree}.
#' 
#' @references Ford, L. R.; Fulkerson, D. R. (1956). "Maximal flow through a 
#' network". Canadian Journal of Mathematics 8: 399.
#' 
#' @examples
#' # Graph
#' nodes <- 1:6
#' arcs <- matrix(c(1,2,1, 1,3,7, 2,3,1, 2,4,3, 2,5,2, 3,5,4, 4,5,1, 4,6,6, 
#'                 5,6,2), byrow = TRUE, ncol = 3)
#' # Maximum flow with Ford-Fulkerson algorithm
#' maxFlowFordFulkerson(nodes, arcs, source.node = 2, sink.node = 6)
#' 
#' @export

maxFlowFordFulkerson <- function(nodes, arcs, directed = FALSE,
                                 source.node = 1, 
                                 sink.node = nodes[length(nodes)]) {
  
  # In case of a undirected graph duplicate and order arcs
  if (!directed) {
    arcs <- rbind(arcs, matrix(c(arcs[, 2], arcs[, 1], arcs[, 3]), ncol = 3))
    arcs <- arcs[order(arcs[, 1], arcs[, 2]), ]
  }
  
  # Previously add a column with flow equal to the capacity of each arc
  if (ncol(arcs) < 4) {
    arcs <- cbind(arcs, arcs[, 3])
  }
  
  # Step 1: found one path between source and sink
  p.arcs <- searchFlowPath(nodes, arcs, source.node, sink.node)$path.arcs
  
  # Step 2: Decrease capacities
  # Repeat until there isn't path between source and sink
  while (nrow(p.arcs) > 0) {
    
    # Decrease capacities to path arcs
    min.flow <- min(p.arcs[, 4])  # minimum capacity
    p.arcs[, 4] <- p.arcs[, 4] - min.flow  # subtract capacities
    
    # Update data in original arc list
    for (i in 1:nrow(p.arcs)) {
      iArcs <- which(arcs[, 1] == p.arcs[i, 1] & arcs[, 2] == p.arcs[i, 2])
      arcs[iArcs, 4] <- p.arcs[i, 4]  # update capacity
    }
    
    # Repeat the search of a path in updated list
    p.arcs <- searchFlowPath(nodes, arcs, source.node, sink.node)$path.arcs
        
  }
  
  
  # Step 3: Determine minimum cut
  stCut <- findstCut(nodes, arcs, source.node, sink.node)
  
  # Step 4: Get maximum flow from the minimum cut
  # Arcs leaving s cut nodes and reaching t cut nodes
  cut.set <- which(arcs[, 1] %in% stCut$s.cut & arcs[, 2] %in% stCut$t.cut)
  # Weight of the cut set equal to maximum flow of the network
  cut.weight <- sum(arcs[cut.set, 3])
  max.flow <- cut.weight
  
  # Build output
  output <- list("s.cut" = stCut$s.cut, "t.cut" = stCut$t.cut,
                 "max.flow" = max.flow)
  return(output)
  
}
#-----------------------------------------------------------------------------#