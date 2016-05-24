#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum spanning trees                                       #
#-----------------------------------------------------------------------------#

# mstEROKruskal ---------------------------------------------------------------
#' ERO rule for minimum cost spanning tree problems with Kruskal's algorithm
#'
#' Given a graph with at least one minimum cost spanning tree, the
#' \code{mstEROKruskal} function divides the cost of the tree among the agents
#' according to the ERO rule.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number from \eqn{1} until the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{mstEROKruskal} returns a matrix with the agents and their
#'  costs.
#' 
#' @references V. Feltkamp, S. H. Tijs, S. Muto, "On the irreducible core and
#' the equal remaining obligation rule of minimum cost extension problems",
#' Mimeo, Tilburg University, 1994.
#' 
#' @seealso The more general function \link{mstRules}.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)
#' # ERO with Kruskal
#' mstEROKruskal(nodes, arcs)

mstEROKruskal <- function(nodes, arcs) {
  
  # Get minimum cost spanning tree and save it
  mst <- getMinimumSpanningTree(nodes, arcs, algorithm = "Kruskal", 
                                show.graph = FALSE, show.data = FALSE)
  msTree <- mst$tree.arcs
  
  # Subtract one to treat nodes as agents
  msTree[, 1:2] <- msTree[, 1:2] - 1
  # Remove column names from tree
  colnames(msTree) <- NULL
  
  # Agents
  players <- unique(c(msTree[, 1], msTree[, 2]))
  # Remove source node
  players <- players[-which(players == 0)]
  # Order players
  players <- players[order(players)]
  # Number of players
  n <- length(players)
  
  # Obligations matrices related to each arc
  obBefMat <- matrix(1, ncol = n)  # previous obligations (all start with one)
  obAftMat <- matrix(ncol = n)[-1, ]  # later obligations
  obArcMat <- matrix(ncol = n)[-1, ]  # obligations to pay for each arch
  
  # Logic value to specify if source are with agents
  sourcePlayers <- FALSE
  
  # Rebuild tree arc by arc
  arcsTree <- c()
  
  # Iterate in each arc of the tree
  for (i in 1:nrow(msTree)) {
    
    # Add arc i
    arcsTree <- rbind(arcsTree, msTree[i, ])
    
    # Affected agents by the arc
    arcPlayers <- unique(c(arcsTree[i, 1], arcsTree[i, 2]))
    if (0 %in% arcPlayers) {
      # Remove source and mark arc
      arcPlayers <- arcPlayers[-which(arcPlayers == 0)]
      sourcePlayers <- TRUE
    }
    
    # Check if there is more affected agents reviewing previews connections
    k <- 1  # start by first agent
    while (k <= length(arcPlayers)) {
      
      if (any(arcPlayers[k] %in% matrix(arcsTree[-i, 1:2], ncol = 2))) {
        # In which arcs appear agent k
        arcs <- which(arcsTree[-i, 1] %in% arcPlayers[k] 
                      | arcsTree[-i, 2] %in% arcPlayers[k])
        # Extract new affected agent
        newPlayers <- unique(c(arcsTree[arcs, 1], arcsTree[arcs, 2]))
        # Add new agent to affected
        arcPlayers <- unique(c(arcPlayers, newPlayers))
        # Remove source and mark arc
        if (0 %in% arcPlayers) {
          arcPlayers <- arcPlayers[-which(arcPlayers == 0)]
          sourcePlayers <- TRUE
        }
      }
      # Continue checking next agent
      k <- k + 1
      length(arcPlayers)
    }    
    
    # Obligations to distribute
    if (sourcePlayers) {
      # 0 if the agents connects to the source
      obligations <- 0
      sourcePlayers <- FALSE
    } else {
      # Otherwise divides between agents
      obligations <- 1/length(arcPlayers)
    }
    
    # Previous obligations
    obBefMat[i, ]
    # Next obligations (modify obligations of affected agents)
    obAftMat <- rbind(obAftMat, obBefMat[i, ])  # start with previous
    obAftMat[i, arcPlayers] <- obligations  # modify
    # Compute how many have to pay each agent
    obArcMat <- rbind(obArcMat, obBefMat[i, ] - obAftMat[i, ])
    
    # Save new previous obligations
    obBefMat <- rbind(obBefMat, obAftMat[i, ])
    
    # Repeat process with next arc of the tree
    
  }
  
  # Computes cost of each arc for every agent
  costs <- c()
  for (i in 1:n) {
    costs <- c(costs, sum(msTree[, 3] * obArcMat[, i]))
  }
  
  # Matrix to store the corresponding costs for each agent
  EROMat <- matrix(c(players + 1, costs), ncol = 2)
  
  # Set appropriate column names
  colnames(EROMat) <- c("agent", "cost")
  
  # Return cost matrix
  return(EROMat)
  
}
#-----------------------------------------------------------------------------#