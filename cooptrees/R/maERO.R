#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum cost arborescences problems                          #
#-----------------------------------------------------------------------------#

 # maERO ----------------------------------------------------------------------
#' ERO rule for minimum cost arborescence problems
#'
#' Given a graph with a minimum cost arborescence, the \code{maERO} function
#' divides the cost of the arborescence among the agents according to the ERO
#' rule. For that purpose, the irreducible form of the problem is obtained. The
#' ERO rule is just the Shapley value of the cooperative game associated with
#' the irreducible form.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{maERO} returns a matrix with the agents and their costs.
#' 
#' @references B. Dutta and D. Mishra, "Minimum cost arborescences", Games and
 #' Economic Behavior, vol. 74, pp. 120-143, Jan. 2012.
#' 
#' @seealso The more general function \link{maRules}.
#' 
#' @examples
#' # Graphs
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,7, 1,3,6, 1,4,4, 2,3,8, 2,4,6, 3,2,6,
#'                  3,4,5, 4,2,5, 4,3,7), ncol = 3, byrow = TRUE)
#' # ERO       
#' maERO(nodes, arcs)

maERO <- function(nodes, arcs) {
  
  # Obtain irreducible form
  irr.form <- maIrreducible(nodes, arcs)
  
  # Cooperative game associated with the irreducible form
  coop.game <- maGames(nodes, irr.form, show.data = FALSE)$values
  
  # Agents
  players <- nodes[-1]
  
  # Compute Shaple value of the game
  costs <- shapleyValue(n = length(players), v = coop.game)$value
  
  # Matrix to store the corresponding costs for each agent
  EROMat <- matrix(c(players, costs), ncol = 2)
  
  # Set appropriate column names
  colnames(EROMat) <- c("agent", "cost")
  
  # Return cost matrix
  return(EROMat)
  
}
#-----------------------------------------------------------------------------#