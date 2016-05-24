#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Shortest Path Tree Problems                                                 #
#-----------------------------------------------------------------------------#

# spTreeDijkstra --------------------------------------------------------------
#' Shortest path tree with Dijkstra's algorithm
#'
#' The \code{spTreeDijkstra} function computes the shortest path tree of an
#' undirected or directed graph with Dijkstra's algorithm.
#' 
#' @details Dijkstra's algorithm was developed by the computer scientist Edsger
#' Dijkstra in 1956 and published in 1959. This is an algorithm that can 
#' computes a shortest path tree from a given source node to the others nodes
#' that make a connected graph, directed or not, with non-negative weights.
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
#' @return \code{spTreeDijkstra} returns a list with:
#' tree.nodes vector containing the nodes of the shortest path tree.
#' tree.arcs matrix containing the list of arcs of the shortest path tree.
#' stages number of stages required.
#' distances vector with distances from source to other nodes
#' 
#' @references Dijkstra, E. W. (1959). "A note on two problems in connexion 
#' with graphs". Numerische Mathematik 1, 269-271.
#' 
#' @seealso A more general function \link{getShortestPathTree}.

spTreeDijkstra <- function(nodes, arcs, source.node = 1, directed = TRUE) {
  
  # Only works with non-negative weights
  if (any(arcs[, 3] < 0)) {
    stop("Dijkstra can't compute negative weights")
  }
  
  # In case of a undirected graph duplicate arcs
  if (directed == FALSE) {
    arcs <- rbind(arcs, matrix(c(arcs[, 2], arcs[, 1], arcs[, 3]), ncol = 3))
  }
  
  # Initialize
  tree.nodes <- nodes[source.node]  # nodes already in the tree
  w <- rep(0, length(nodes))  # distances from source to each node
  tree.arcs <- matrix(ncol = 3)[-1, ]  # matrix to save shortest path tree
  
  stages <- 0  # initialize counter
  
  # Iterate until all the nodes are in the tree
  while (length(tree.nodes) < length(nodes)) {
    
    # Select arcs between checked and unchecked nodes
    k <- which(arcs[, 1] %in% tree.nodes & 
                 arcs[, 2] %in% nodes[-which(nodes %in% tree.nodes)])
    kArcs <- matrix(arcs[k, ], ncol = 3)
    
    # Get minimum weight and one arc with it
    minW <- min(kArcs[, 3] + w[kArcs[, 1]])
    min.arc <- matrix(kArcs[which(kArcs[, 3] + 
                                   w[kArcs[, 1]] == minW), ], ncol = 3)
    
    # Save selected arc and distance
    tree.arcs <- rbind(tree.arcs, min.arc[1, ])
    tree.nodes <- c(tree.nodes, min.arc[1, 2])
    w[min.arc[1, 2]] <- w[min.arc[1, 2]] + minW
    
    stages <- stages + 1  # counter
  }
  
  # Column names
  colnames(tree.arcs) <- c("head", "tail", "weight")
  
  output <- list("tree.nodes" = tree.nodes, "tree.arcs" = tree.arcs,
                 "stages" = stages, "distances" = w)
  return(output)
  
}
#-----------------------------------------------------------------------------#