#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum cost arborescences problems                          #
#-----------------------------------------------------------------------------#

# maBird ----------------------------------------------------------------------
#' Bird rule for minimum cost arborescence problems
#'
#' Given a graph with a minimum cost arborescence, the \code{maBird} function
#' divides the cost of this arborescence among the agents. For that purpose it,
#' uses the Bird rule, where each agent pays the cost of the last arc that 
#' connects him to the source.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{maBird} returns a matrix with the agents and their costs.
#' 
#' @references B. Dutta and D. Mishra, "Minimum cost arborescences", Games and
#' Economic Behavior, vol. 74, pp. 120-143, Jan. 2012.
#' 
#' @seealso The more general function \link{maRules}.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,7, 1,3,6, 1,4,4, 2,3,8, 2,4,6, 3,2,6,
#'                  3,4,5, 4,2,5, 4,3,7), ncol = 3, byrow = TRUE)
#' # Bird                  
#' maBird(nodes, arcs)

maBird <- function(nodes, arcs) {
  
  # Get minimum cost arborescence and save it
  arbor <- getMinimumArborescence(nodes, arcs, show.graph = FALSE,
                                  show.data = FALSE)
  arborescence <- arbor$tree.arcs
  
  # Each agent pay the cost of the arc that connects it to the arborescence
  agents <- arborescence[, 2]
  costs <- arborescence[, 3]
  
  # Matrix to store the corresponding costs for each agent
  birdMat <- matrix(c(agents, costs), ncol = 2)
  birdMat <- matrix(birdMat[order(birdMat[, 1]), ], ncol = 2)  # reorder by agents
  
  # Set appropriate column names
  colnames(birdMat) <- c("agent", "cost")
  
  # Return cost matrix
  return(birdMat)
  
}
#-----------------------------------------------------------------------------#