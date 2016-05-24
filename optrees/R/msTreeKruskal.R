#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Spanning Tree Problem                                               #
#-----------------------------------------------------------------------------#

# msTreeKruskal ---------------------------------------------------------------
#' Minimum cost spanning tree with Kruskal's algorithm
#'
#' \code{msTreeKruskal} computes a minimum cost spanning tree of an undirected 
#' graph with Kruskal's algorithm.
#' 
#' @details Kruskal's algorithm was published for first time in 1956 by 
#' mathematician Joseph Kruskal. This is a greedy algorithm that finds a 
#' minimum cost spanning tree in a connected weighted undirected graph by
#' adding, without forming cycles, the minimum weight arc of the graph at each
#' stage.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' 
#' @return \code{msTreeKruskal} returns a list with:
#' tree.nodes vector containing the nodes of the minimum cost spanning tree.
#' tree.arcs matrix containing the list of arcs of the minimum cost spanning tree.
#' stages number of stages required.
#' stages.arcs stages in which each arc was added.
#' 
#' @references Kruskal, Joshep B. (1956), "On the Shortest Spanning Subtree of a
#' Graph and the Traveling Salesman Problem", Proceedings of the American 
#' Mathematical Society, Vol. 7, No. 1 (Feb., 1956), pp. 48-50
#' 
#' @seealso A more general function \link{getMinimumSpanningTree}.
#' 
#' @export

msTreeKruskal <- function(nodes, arcs) {
  
  # Order arcs by weight
  arcs <- matrix(arcs[order(arcs[, 3]), ], ncol = 3)
  
  # Components
  components <- matrix(c(nodes, nodes), ncol = 2)
  
  # Initialize tree with first arc
  tree.arcs <- matrix(ncol = 3)[-1, ]
  
  stages <- 0  # initialize counter
  stages.arcs <- c()  # vector to store stage number in wich each arc was added
  
  # Start with first arc
  i <- 1
  # Repeat until we have |N|-1 arcs
  while(nrow(tree.arcs) < length(nodes) - 1) {
    
    # Select arc
    min.arc <- arcs[i, ]
    
    # Check components of the two nodes of selected arc
    iComp <- components[components[, 1] == min.arc[1], 2]
    jComp <- components[components[, 1] == min.arc[2], 2]
    if (iComp != jComp) {
      # Add arc to msTree
      tree.arcs <- rbind(tree.arcs, min.arc)
      # Merge components
      components[components[, 2] == jComp, 2] <- iComp
    }
    
    stages <- stages + 1  # counter
    # Save in which stage an arc was added to the tree and update
    stages.arcs <- c(stages.arcs,
                     rep(stages, nrow(tree.arcs) - length(stages.arcs)))
    # Continue with next arc
    i <- i + 1
    
  }
  
  # Column names
  colnames(tree.arcs) <- c("ept1", "ept2", "weight")
  # Remove row names
  rownames(tree.arcs) <- NULL
  
  output <- list("tree.nodes" = nodes, "tree.arcs" = tree.arcs,
                 "stages" = stages, "stages.arcs" = stages.arcs)
  return(output)
  
}
#-----------------------------------------------------------------------------#