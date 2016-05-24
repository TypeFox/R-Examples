#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum cost arborescences problems                          #
#-----------------------------------------------------------------------------#

# maPessimistic ---------------------------------------------------------------
#' Pessimistic game associated with minimum cost arborescences
#' 
#' Given a graph with at least one minimum cost arborescence, the
#' \code{maPessimistic} function builds the pessimistic game.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{maPessimistic} returns a vector with the characteristic 
#' function of the pessimitic game.
#' 
#' @references B. Dutta and D. Mishra, "Minimum cost arborescences", Games and
#' Economic Behavior, vol. 74, pp. 120-143, Jan. 2012.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,7, 1,3,6, 1,4,4, 2,3,8, 2,4,6, 3,2,6,
#'                  3,4,5, 4,2,5, 4,3,7), ncol = 3, byrow = TRUE)
#' # Pessimistic game
#' maPessimistic(nodes, arcs)

maPessimistic <- function(nodes, arcs) {
  
  # Initialize
  vS <- c()  # values to save
  players <- c(1:(length(nodes)-1))  # agents
  coalitions <- c()  # vector to store
  
  # Work with cost matrix
  Cmat <- ArcList2Cmat(nodes, arcs, directed = TRUE)
  
  # Iterate for coalitions of diferent size
  for (i in 1:max(players)) {
    
    # If all agents are on the coalition solve directly
    if (i == length(players)) {
      
      # Get minimum cost arborescence
      temparbor <- getMinimumArborescence(nodes, arcs, show.graph = FALSE,
                                          show.data = FALSE)
      
      # Check if there is an arborescence
      if (class(temparbor) == "logical") {
        # If not the cost is 0
        vS <- c(vS, 0)
      } else {
        # Otherwise compute cost of arborescence and save it
        vS <- c(vS, sum(temparbor$tree.arcs[, 3]))
      }
      # Save coalition
      coalitions <- c(coalitions, paste(players, sep="", collapse=","))
      
    } else {
      
      # Check every possible combination of i agents
      S <- combn(length(players), i)
      
      # Iterate for each coalitions with size i
      for (j in 1:ncol(S)) {
        
        # Get coalition to check
        Sactual <- S[, j]
        coalitions <- c(coalitions, paste(Sactual, sep="", collapse=","))
        # Select nodes and arcs from the actual coalition
        Snodes <- c(0, Sactual) + 1  # coalition nodes + source
        # Cost matrix of the coalition + source
        sCmat <- Cmat[-which(!nodes %in% Snodes), -which(!nodes %in% Snodes)]
        tempnodes <- 1:nrow(sCmat)
        temparcs <- Cmat2ArcList(tempnodes, sCmat)
        temparbor <- getMinimumArborescence(tempnodes, temparcs, 
                                            show.graph = FALSE,
                                            show.data = FALSE)
        
        # Check if there is an arborescence
        if (class(temparbor)=="logical") {
          # If not the cost is 0
          vS <- c(vS, 0)
        } else {
          # Otherwise compute cost of arborescence and save it
          cost <- sum(temparbor$tree.arcs[,3])
        }  
        
        # Save costs
        vS <- c(vS, cost)
        
      }
    }
  }
  
  # Returns vectors with coalitions and the characteristic function
  return(list(coalitions = coalitions, values = vS))
  
}
#-----------------------------------------------------------------------------#