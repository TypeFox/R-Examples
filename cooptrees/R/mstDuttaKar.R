#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum spanning trees                                       #
#-----------------------------------------------------------------------------#

# mstDuttaKar -----------------------------------------------------------------
#' Dutta-Kar rule for minimum cost spanning tree problems
#' 
#' Given a graph with at least one minimum cost spanning tree, the
#' \code{mstDuttaKar} function divides the cost of the tree obtained with
#' Prim's algorithm among the agents according to the Dutta-Kar rule. This rule
#' specifies that each agent chooses to pay the minimum cost between the last
#' arc that connects him to the source and the cost that rejects his successor.
#' The order is set by Prim's algorithm.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{mstDuttaKar} returns a matrix with the agents and their costs.
#' 
#' @references B. Dutta and A. Kar, "Cost monotonicity, consistency and minimum
#' cost spanning tree games", Games and Economic Behavior, vol. 48,
#' pp. 223-248, Aug. 2004.
#' 
#' @seealso The more general function \link{mstRules}.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)
#' # Dutta-Kar
#' mstDuttaKar(nodes, arcs)

mstDuttaKar <- function(nodes, arcs) {
  
  # Get minimum spanning tree and save it
  mst <- getMinimumSpanningTree(nodes, arcs, algorithm = "Prim", 
                                show.graph = FALSE, show.data = FALSE)
  msTree <- mst$tree.arcs
  
  # Agents in the tree
  players <- c()
  for (i in 1:nrow(msTree)) {
    players <- c(players, msTree[i, 1], msTree[i, 2])
  }
  players <- unique(players)
  # Remove source node and subtract one to operate with indexs in R
  players <- players[-which(players == 1)]
  # Vector to store the costs
  costs <- c()
  
  # Initialize with cost of first arc
  costPend <- msTree[1, 3]; costPend
  
  # Each agent choose between pending cost and the cost of de arc leaving it
  for (i in 1:length(players)) {
    
    # Iterate among every agent
    arcPlayer <- i + 1; arcPlayer  # arc leaving i
    
    if (arcPlayer > nrow(msTree)) {
      # If there is no arc leaving the i agent assumes pending cost
      costPlayer <- costPend
      # New costPend of the next arc
    } else {
      # If there is arcs leaving the i agent assumes minimum cost
      costPlayer <- min(c(msTree[arcPlayer, 3], costPend)); costPlayer  # minimum cost
      costPend <- max(msTree[arcPlayer, 3], costPend); costPend  # new pending cost
    }
    
    # Assign cost to agent
    costs <- c(costs, costPlayer)
  }
  
  # Matrix to store the corresponding costs for each agent
  duttakarMat <- matrix(c(players, costs), ncol = 2)
  
  # Set appropriate column names
  colnames(duttakarMat) <- c("agent", "cost")
  
  # Return cost matrix
  return(duttakarMat)
  
}
#-----------------------------------------------------------------------------#