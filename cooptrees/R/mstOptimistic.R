#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum spanning trees                                       #
#-----------------------------------------------------------------------------#

# mstOptimistic ---------------------------------------------------------------
#' Optimistic game of a minimum cost spanning tree problem
#'
#' Given a graph with at least one minimum cost spanning tree, the
#' \code{mstOptimistic} function builds the optimistic game associated with it.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{mstOptimistic} returns a vector with the characteristic
#' function of the optimistic game.
#' 
#' @references G. Bergantiños and J. J. Vidal-Puga, "The optimistic TU game in
#' minimum cost spanning tree problems", International Journal of Game Theory,
#' vol. 36, pp. 223-239, Feb. 2007.
#' 
#' @seealso The more general function \link{mstGames}.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)
#' # Optimistic game
#' mstOptimistic(nodes, arcs)

mstOptimistic <- function(nodes, arcs) {
  
  # Initialize
  vS <- c()  # vector to store characteristic function
  players <- nodes[-1] - 1  # agents without source node
  coalitions <- c()  # coalitions
  
  # Iterate for coalitions of diferent size
  for (i in 1:max(players)) {
    
    # If all agents are on the coalition solve directly
    if (i == length(players)) {
      
      # Get minimum spanning tree
      mst <- getMinimumSpanningTree(nodes, arcs, algorithm = "Prim",
                                    show.graph = FALSE, show.data = FALSE)
      # Save cost and coalition
      vS <- c(vS, mst$weight)
      coalitions <- c(coalitions, paste(players, sep="", collapse=","))
      
    } else {  
      
      # Check every possible combination of i agents
      S <- combn(length(players), i)
      
      # Iterate for each coalitions with size i
      for (j in 1:ncol(S)) {
        
        # Get coalition to check
        Sactual <- S[, j]
        # Select nodes and arcs from de actual coalition
        Snodes <- c(0, Sactual) + 1  # coalition nodes + source
        Sarcs <- matrix(arcs[which(arcs[, 1] %in% Snodes
                                   & arcs[, 2] %in% Snodes), ], ncol = 3)
        
        # Recalculate costs between coalition agents and the source
        sourceArcs <- which(Sarcs[, 1] == 1 | Sarcs[, 2] == 1)  # arcs
        arcsPl <- unique(c(Sarcs[sourceArcs, 1:2]))  # affected agents
        arcsPl <- arcsPl[-which(arcsPl == 1)]  # without the source
        
        # Check arcs that leaving from each coalition agents to others
        for (k in 1:length(arcsPl)) {
          
          # Arcs leaving or entering in agent k
          kArcs <- which(arcs[, 1] == arcsPl[k] | arcs[, 2] == arcsPl[k])
          # Select only arcs leaving or entering agents outside the coalition
          extPlayers <- players[-which(players %in% (arcsPl - 1))] + 1
          kextArcs <- which(arcs[kArcs, 1] %in% c(1, extPlayers)
                            | arcs[kArcs, 2] %in% c(1, extPlayers))
          extArcs <- matrix(arcs[kArcs, ][kextArcs, ], ncol = 3)
          # Get minimum cost from that arcs
          minCost <- min(extArcs[, 3])
          # Add to arc that links agent k with the source
          Sarcs[which(Sarcs[, 1] == 1 
                      & Sarcs[, 2] == arcsPl[k]), 3] <- minCost
          
        }
        
        validGraph <- checkGraph(Snodes, Sarcs)
        if (validGraph) {
          # Get minimum spanning tree
          mst <- getMinimumSpanningTree(Snodes, Sarcs, algorithm = "Prim",
                                        show.graph = FALSE, show.data = FALSE)
          vWeight <- mst$weight
        } else {
          vWeight <- NA
        }
        # Save cost and coalition
        vS <- c(vS, vWeight)
        coalitions <- c(coalitions, paste(Sactual, sep="", collapse=","))
      }
    }
  }
  
  # Returns vectors with coalitions and the characteristic function
  return(list(coalitions = coalitions, values = vS))
  
}
#-----------------------------------------------------------------------------#