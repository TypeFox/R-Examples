#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum spanning trees                                       #
#-----------------------------------------------------------------------------#

# mstKar ----------------------------------------------------------------------
#' Kar rule for minimum cost spanning tree problems
#' 
#' Given a graph with at least one minimum cost spanning tree, the
#' \code{mstKar} function divides the cost of the tree among the agents
#' according to the Kar rule. That rule is obtained with the Shapley value of
#' the pessimistic game.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{mstKar} returns a matrix with the agents and their costs.
#' 
#' @references A. Kar, "Axiomatization of the Shapley Value on Minimum Cost
#' Spanning Tree Games", Games and Economic Behavior, vol. 38, pp. 265-277,
#' Feb. 2002.
#' 
#' @seealso The more general function \link{mstRules}.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)
#' # Kar
#' mstKar(nodes, arcs)

mstKar <- function(nodes, arcs) {
  
  # Get pessimistic game
  pesGame <- mstPessimistic(nodes, arcs)
  
  # Agents in the graph
  players <- nodes[-1]
  
  # Values and list of coalitions of the pessimistic game
  S <- pesGame$coalitions
  v <- pesGame$values
  
  # Compute Shapley value
  costs <- shapleyValue(n = length(players), S = S, v = v)$value
  
  # Matrix to store the corresponding costs for each agent
  karMat <- matrix(c(players, costs), ncol = 2)
  
  # Set appropriate column names
  colnames(karMat) <- c("agent", "cost")
  
  # Return cost matrix
  return(karMat)
  
}
#-----------------------------------------------------------------------------#