#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum cost arborescences                                   #
#-----------------------------------------------------------------------------#

# maGames ---------------------------------------------------------------------
#' Cooperative game associated with minimum cost arborescences
#' 
#' Given a graph with at least one minimum cost arborescence the \code{maGames}
#' function builds the cooperative game associated with it.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param game denotes the game to be obtained, known as the "pessimistic"
#' game.
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). By default its value is
#' \code{TRUE}.
#'
#' @return \code{maGames} returns a vector with the characteristic fuction of
#' the cooperative game associated with the graph.
#' 
#' @references B. Dutta and D. Mishra, "Minimum cost arborescences", Games and
#' Economic Behavior, vol. 74, pp. 120-143, Jan. 2012.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,7, 1,3,6, 1,4,4, 2,3,8, 2,4,6, 3,2,6,
#'                  3,4,5, 4,2,5, 4,3,7), ncol = 3, byrow = TRUE)
#' # Cooperative games
#' maGames(nodes, arcs, game = "pessimistic")

maGames <- function(nodes, arcs, game = "pessimistic", show.data = TRUE) {
  
  # Obtain game according to selection
  if (game == "pessimistic") {
    coopGames <- maPessimistic(nodes, arcs)
  } else {
    warning("Unknown game. Choose pessimistic.")
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