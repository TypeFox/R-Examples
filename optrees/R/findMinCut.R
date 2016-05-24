#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimum Cut Tree Problems                                                   #
#-----------------------------------------------------------------------------#

# findMinCut ------------------------------------------------------------------
#' Finds the minimum cut of a given graph
#' 
#' The \code{findMinCut} function can find the minimum cut of a given graph.
#' For that, this function computes the maximum flow of the network and applies
#' the max-flow min-cut theorem to determine the cut with minimum weight 
#' between the source and the sink nodes.
#' 
#' @details The max-flow min-cut theorem proves that, in a flow network, the
#' maximum flow between the source node and the sink node and the weight of any
#' minimum cut between them is equal.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param algorithm denotes the algorithm used to compute the maximum flow:
#' "Ford-Fulkerson".
#' @param source.node number pointing to the source node of the graph. It's node
#' \eqn{1} by default.
#' @param sink.node number pointing to the sink node of the graph. It's the
#' last node by default.
#' @param directed logical value indicating whether the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#' 
#' @return \code{findMinCut} returns a list with:
#' s.cut vector with the nodes of the \code{s} cut.
#' t.cut vector with the nodes of the \code{t} cut.
#' maxFlow value with the maximum flow in the flow network.
#' cut.set list of arcs of the cut set founded.
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{ghTreeGusfield} and \link{getMinimumCutTree}.
#' 
#' @examples
#' # Graph
#' nodes <- 1:6
#' arcs <- matrix(c(1,2,1, 1,3,7, 2,3,1, 2,4,3, 2,5,2, 3,5,4, 4,5,1, 4,6,6, 
#'                 5,6,2), byrow = TRUE, ncol = 3)
#' # Find minimum cut
#' findMinCut(nodes, arcs, source.node = 2, sink.node = 6)
#' 
#' @export

findMinCut <- function(nodes, arcs, algorithm = "Ford-Fulkerson", 
                       source.node = 1, sink.node = nodes[length(nodes)],
                       directed = FALSE) {
  
  # In case of a undirected graph duplicate and order arcs
  if (!directed) {
    arcs <- rbind(arcs, matrix(c(arcs[, 2], arcs[, 1], arcs[, 3]), ncol = 3))
    arcs <- arcs[order(arcs[, 1], arcs[, 2]), ]
  }
  
  # Add a column with flow equal to the capacity of each arc
  arcs <- cbind(arcs, arcs[, 3])
  colnames(arcs) <- NULL
  
  # Get minimum cut with selected algorithm
  if (algorithm == "Ford-Fulkerson") {
    min.cut <- maxFlowFordFulkerson(nodes, arcs, TRUE, source.node, sink.node)
  } else {
    stop("Unknown algorithm")
  }
  
  # Get cut set
  cut.arcs <- which(arcs[, 1] %in% min.cut$s.cut & arcs[, 2] %in% min.cut$t.cut)
  cut.set <- matrix(arcs[cut.arcs, -4], ncol = 3)
  colnames(cut.set) <- c("ept1", "ept2", "weight")
  min.cut <- append(min.cut, list("cut.set" = cut.set))
  
  # Return minimum cut
  return(min.cut)
  
}
#-----------------------------------------------------------------------------#