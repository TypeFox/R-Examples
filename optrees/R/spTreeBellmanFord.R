#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Shortest Path Tree Problems                                                 #
#-----------------------------------------------------------------------------#

# spTreeBellmanFord -----------------------------------------------------------
#' Shortest path tree with Bellman-Ford algorithm
#'
#' The \code{spTreeBellmanFord} function computes the shortest path tree of an
#' undirected or directed graph with the Bellman-Ford algorithm.
#' 
#' @details The Bellman-Ford algorithm gets its name for two of the developers,
#' Richard Bellman y Lester Ford Jr., that published it in 1958 and 1956
#' respectively. The same algorithm also was published independently in 1957
#' by Edward F. Moore.
#'  
#' The Bellman-Ford algorithm can compute the shortest path from a source node
#' to the rest of nodes that make a connected graph, directed or not, with
#' weights that can be negatives. If the graph is connected and there isn't
#' negative cycles, the algorithm always finds a shortest path tree.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param source.node number pointing the source node of the graph. It's node
#' \eqn{1} by default.
#' @param directed logical value indicating whether the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#'
#' @return \code{spTreeBellmanFord} returns a list with:
#' tree.nodes vector containing the nodes of the shortest path tree.
#' tree.arcs matrix containing the list of arcs of the shortest path tree.
#' stages number of stages required.
#' distances vector with distances from source to other nodes
#' 
#' @references Bellman, Richard (1958). "On a routing problem". Quarterly of 
#' Applied Mathematics 16, 87-90.
#' 
#' Ford Jr., Lester R. (1956). Network Flow Theory. Paper P-923. Santa Monica, 
#' California: RAND Corporation.
#' 
#' Moore, Edward F. (1959). "The shortest path through a maze". Proc. Internat.
#' Sympos. Switching Theory 1957, Part II. Cambridge, Mass.: Harvard Univ.
#' Press. pp. 285-292.
#' 
#' @seealso A more general function \link{getShortestPathTree}.

spTreeBellmanFord <- function(nodes, arcs, source.node = 1, directed = TRUE) {
  
  # In case of a undirected graph duplicate arcs
  if (directed == FALSE) {
    arcs <- rbind(arcs, matrix(c(arcs[, 2], arcs[, 1], arcs[, 3]), ncol = 3))
  }
  
  # Initialize
  n <- length(nodes)  # order of the graph
  # Initialize distances
  tree.nodes <- matrix(c(nodes, rep(Inf, n), numeric(n)), ncol = 3)
  tree.nodes[source.node, 2] <- 0
  tree.arcs <- matrix(ncol = 3)[-1, ] # matrix to save shortest path tree
  
  stages <- 0  # initialize counter
  
  # Relax arcs repeatdly
  for (i in 1:(n-1)) {
    for (i in 1:nrow(arcs)) {
      if (tree.nodes[arcs[i, 2], 2] > tree.nodes[arcs[i, 1], 2] + arcs[i, 3]) {
        # If distance stored in reaching node is greater than distance from
        # precedent node plus the weight of the arc, save new distance
        tree.nodes[arcs[i, 2], 2] <- tree.nodes[arcs[i, 1], 2] + arcs[i, 3]
        tree.nodes[arcs[i, 2], 3] <- arcs[i, 1]
      }      
    }
    stages <- stages + 1
  }
  
  # Check for negative-weight cycles
  for (i in 1:nrow(arcs)) {
    if (tree.nodes[arcs[i, 2], 2] > tree.nodes[arcs[i, 1], 2] + arcs[i, 3]) {
      stop("There is a negative-weight cycle")
    }
  }
  
  colnames(tree.nodes) <- c("node", "dist", "pred")
  
  # Rebuilt shortest path tree
  for (i in 1:nrow(tree.nodes)) {
    k <- which(arcs[, 1] == tree.nodes[i, 3] & arcs[, 2] == tree.nodes[i, 1])
    spArc <- arcs[k, ]
    tree.arcs <- rbind(tree.arcs, spArc)
  }
  
  # Reorder arcs
  tree.arcs <- tree.arcs[order(tree.arcs[, 1], tree.arcs[, 2]), ]
  rownames(tree.arcs) <- NULL
  # Column names
  colnames(tree.arcs) <- c("head", "tail", "weight")
  
  output <- list("tree.nodes" = nodes, "tree.arcs" = tree.arcs,
                 "stages" = stages, "distances" = tree.nodes[, 2])
  return(output)
  
}
#-----------------------------------------------------------------------------#