#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Spanning Tree Problem                                               #
#-----------------------------------------------------------------------------#

# msTreePrim ------------------------------------------------------------------
#' Minimum cost spanning tree with Prim's algorithm
#'
#' \code{msTreePrim} computes a minimum cost spanning tree of an undirected 
#' graph with Prim's algorithm.
#' 
#' @details Prim's algorithm was developed in 1930 by the mathematician Vojtech 
#' Jarnik, later proposed by the computer scientist Robert C. Prim
#' in 1957 and rediscovered by Edsger Dijkstra in 1959. This is a greedy
#' algorithm that can find a minimum spanning tree in a connected, weighted and
#' undirected graph by adding recursively minimum cost arcs leaving visited
#' nodes.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param start.node number associated with the first node in Prim's algorithm.
#' By default, node \eqn{1} is the first node.
#' 
#' @return \code{msTreePrim} returns a list with:
#' tree.nodes vector containing the nodes of the minimum cost spanning tree.
#' tree.arcs matrix containing the list of arcs of the minimum cost spanning
#' tree.
#' stages number of stages required.
#' stages.arcs stages in which each arc was added.
#' 
#' @references Prim, R. C. (1957), "Shortest Connection Networks And Some
#' Generalizations", Bell System Technical Journal, 36 (1957), pp. 1389-1401
#' 
#' @seealso A more general function \link{getMinimumSpanningTree}.
#' 
#' @export

msTreePrim <- function(nodes, arcs, start.node = 1) {
  
  # Duplicate the arcs
  arcs <- rbind(arcs, matrix(c(arcs[, 2], arcs[, 1], arcs[, 3]), ncol = 3))
  
  # Initialize with empty tree and start.node
  tree.arcs <- matrix(ncol = 3)[-1, ]
  tree.nodes <- nodes[nodes == start.node]
  
  stages <- 0  # initialize counter
  stages.arcs <- c()  # vector to store stage number in wich each arc was added
  
  # Iterate until every node was checked
  while (length(tree.nodes) < length(nodes)) {
    
    # Arcs leaving tree nodes and reaching unadded nodes
    k <- which(arcs[, 1] %in% tree.nodes &
                 arcs[, 2] %in% nodes[-which(nodes %in% tree.nodes)])
    validArcs <- matrix(arcs[k, ], ncol = 3)
    
    # Valid arcs with minimum cost
    l <- which(validArcs[, 3] == min(validArcs[, 3]))
    min.arc <- matrix(validArcs[l, ], ncol = 3)
    
    # Save one arc in minimum tree
    tree.arcs <- rbind(tree.arcs, min.arc[1, ])
    
    # Update tree nodes
    tree.nodes <- c(tree.nodes, min.arc[1, 2])
    
    stages <- stages + 1  # counter
    # Save in which stage an arc was added to the tree and update
    stages.arcs <- c(stages.arcs,
                     rep(stages, nrow(tree.arcs) - length(stages.arcs)))
    
  }
  
  # Column names
  colnames(tree.arcs) <- c("ept1", "ept2", "weight")
  
  output <- list("tree.nodes" = tree.nodes, "tree.arcs" = tree.arcs,
                 "stages" = stages, "stages.arcs" = stages.arcs)
  return(output)
  
}
#-----------------------------------------------------------------------------#