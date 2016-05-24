#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum spanning trees                                       #
#-----------------------------------------------------------------------------#

# mstPessimistic --------------------------------------------------------------
#' Pessimistic game from a minimum cost spanning tree problem
#'
#' Given a graph with at least one minimum cost spanning tree, the
#' \code{mstPessimistic} function builds the pessimistic game.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{mstPessimistic} returns a vector with the characteristic
#' function of the pessimistic game.
#' 
#' @references C. G. Bird, "On Cost Allocation for a Spanning Tree: A Game
#' Theoretic Approach", Networks, no. 6, pp. 335-350, 1976.
#' 
#' @seealso The more general function \link{mstGames}.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)
#' # Pessimistic game
#' mstPessimistic(nodes, arcs)

mstPessimistic <- function(nodes, arcs) {
  
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