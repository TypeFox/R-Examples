#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum spanning trees                                       #
#-----------------------------------------------------------------------------#

# mstIrreducible --------------------------------------------------------------
#' Irreducible form for a minimum cost spanning tree problem
#' 
#' Given a graph with at least one minimum cost spanning tree, the
#' \code{mstIrreducible} function obtains the irreducible form.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{mstIrreducible} returns a matrix with the list of arcs of the
#' irreducible form.
#' 
#' @references C. G. Bird, "On Cost Allocation for a Spanning Tree: A Game
#' Theoretic Approach", Networks, no. 6, pp. 335-350, 1976.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)                 
#' # Irreducible form
#' mstIrreducible(nodes, arcs)

mstIrreducible <- function(nodes, arcs) {
  
  # Get mst
  mst <- getMinimumSpanningTree(nodes, arcs, algorithm = "Prim",
                                show.data = FALSE, show.graph = FALSE)
  msTree <- mst$tree.arcs
  # Duplicate arcs
  msTree <- rbind(msTree, matrix(c(msTree[, 2], msTree[, 1],
                                   msTree[, 3]), ncol = 3))
  
  # Same arcs than original graph
  irreducibleForm <- arcs
  
  # Check arc by arc
  for (i in 1:nrow(arcs)) {
    
    # i and j nodes of the arc i
    iNode <- irreducibleForm[i, 1]
    jNode <- irreducibleForm[i, 2]
    # Get walk from i to j in msTree
    ijWalk <- searchWalk(nodes, msTree,
                         start.node = iNode, end.node = jNode)$walk.arcs
    
    # Maximum cost from the walk
    maxCost <- max(ijWalk[, 3])
    # Change cost to the arc that connects i and j
    irreducibleForm[i, 3] <- maxCost
    
  }
  
  # Column names
  colnames(irreducibleForm) <- c("ept1", "ept2", "weight")
  
  # Return irreducible form
  return(irreducibleForm)
  
}
#-----------------------------------------------------------------------------#