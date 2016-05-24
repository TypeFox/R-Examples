#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum spanning trees                                       #
#-----------------------------------------------------------------------------#

# mstGames --------------------------------------------------------------------
#' Cooperative games from minimum cost spanning tree problems
#' 
#' @description Given a graph with at least one minimum cost spanning tree,
#' \code{mstGames} builds both cooperative games: the pessimistic and the
#' optimistic game.
#' 
#' The pessimistic game associated with a minimum cost spanning
#' tree problem is a cooperative game in which every coalition of agents
#' obtains the minimum cost assuming that the agents outside the coalition
#' are not connected.
#' 
#' The optimistic game associated with with a minimum cost spanning tree
#' problem is a cooperative game in which every coalition of agents obtains
#' the minimum cost assuming that that the agents outside the coalition
#' are connected. Thus, the agents in the coalition can benefit from their
#' connections to the source
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param game denotes the game that we want to obtain: "pessimistic" or 
#' "optimistic".
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). By default its value is
#' \code{TRUE}.
#'
#' @return \code{mstGames} returns a vector with the characteristic fuction of
#' the selected game associated with the graph and prints the result in
#' console.
#' 
#' @references C. G. Bird, "On Cost Allocation for a Spanning Tree: A Game
#' Theoretic Approach", Networks, no. 6, pp. 335-350, 1976.
#' 
#' G. Bergantiños and J. J. Vidal-Puga, "The optimistic TU game in
#' minimum cost spanning tree problems", International Journal of Game Theory,
#' vol. 36, pp. 223-239, Feb. 2007.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)
#' # Cooperative games
#' mstGames(nodes, arcs, game = "pessimistic")
#' mstGames(nodes, arcs, game = "optimistic")

mstGames <- function(nodes, arcs, game, show.data = TRUE) {
  
  # Obtain game according to selection
  if (game == "pessimistic") {
    coopGames <- mstPessimistic(nodes, arcs)
  } else if (game == "optimistic") {
    coopGames <- mstOptimistic(nodes, arcs)
  } else {
    warning("Unknown game. Choose pessimistic or optimisc.")
  }
  
  # Get coalitions and values
  S <- coopGames$coalitions
  v <- coopGames$values
  
  # Console output
  if (show.data) {
    cat("v(S) = ", v, "\n")
  }
  
  # Returns vectors with coalitions and the characteristic function
  return(list(coalitions = S, values = v))
  
}
#-----------------------------------------------------------------------------#