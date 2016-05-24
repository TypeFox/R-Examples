#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum spanning trees                                       #
#-----------------------------------------------------------------------------#

# mstBird ---------------------------------------------------------------------
#' Bird rule for minimum cost spanning tree problems
#'
#' Given a graph with at least one minimum cost spanning tree, the
#' \code{mstBird} function divides the cost of the tree obtained with Prim's
#' algorithm among the agents. For that purpose it uses the Bird rule, where
#' each agent pays the cost of the arc that connects him to the tree source.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{mstBird} returns a matrix with the agents and their costs.
#' 
#' @references C. G. Bird, "On Cost Allocation for a Spanning Tree: A Game
#' Theoretic Approach", Networks, no. 6, pp. 335-350, 1976.
#' 
#' @seealso The more general function \link{mstRules}.
#' 
#' @examples
#' # Graphs
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)
#' # Bird
#' mstBird(nodes, arcs)

mstBird <- function(nodes, arcs) {
  
  # Get minimum spanning tree and save it
  mst <- getMinimumSpanningTree(nodes, arcs, algorithm = "Prim", 
                                show.graph = FALSE, show.data = FALSE)
  msTree <- mst$tree.arcs
  
  # Agents in the tree
  players <- unique(c(msTree[, 1], msTree[, 2]))[-1]
  
  # Each agent pays the cost of the arcs entering in them
  costs <- msTree[, 3]  # so get the arcs in same order
  
  # Matrix to store the corresponding costs for each agent
  birdMat <- matrix(c(players, costs), ncol = 2)
  
  # Set appropriate column names
  colnames(birdMat) <- c("agent", "cost")
  
  # Return cost matrix
  return(birdMat)
  
}
#-----------------------------------------------------------------------------#