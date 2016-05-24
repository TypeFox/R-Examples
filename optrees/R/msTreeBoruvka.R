#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Spanning Tree Problem                                               #
#-----------------------------------------------------------------------------#

# msTreeBoruvka ---------------------------------------------------------------
#' Minimum cost spanning tree with Boruvka's algorithm.
#'
#' \code{msTreeBoruvka} computes a minimum cost spanning tree of an undirected 
#' graph with Boruvka's algorithm.
#' 
#' @details Boruvka's algorithm was firstly published in 1926 by the
#' mathematician Otakar Boruvka. This algorithm works in a connected, weighted 
#' and undirected graph, checking each component and adding the minimum weight
#' arcs that connect the component to other components until one minimum
#' spanning tree is complete.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'  
#' @return \code{msTreeBoruvka} returns a list with:
#' tree.nodes vector containing the nodes of the minimum cost spanning tree.
#' tree.arcs matrix containing the list of arcs of the minimum cost spanning tree.
#' stages number of stages required.
#' stages.arcs stages in which each arc was added.
#' 
#' @references Boruvka, Otakar (1926). "O jistem problemu minimalnim (About a
#' certain minimal problem)". Prace mor. prirodoved. spol. v Brne III (in 
#' Czech, German summary) 3: 37-58.
#' 
#' @seealso A more general function \link{getMinimumSpanningTree}.
#' 
#' @export

msTreeBoruvka <- function(nodes, arcs) {
  
  # Duplicate the arcs
  arcs <- rbind(arcs, matrix(c(arcs[, 2], arcs[, 1], arcs[, 3]), ncol = 3))
  
  # Components
  components <- matrix(c(nodes, nodes), ncol = 2)
  
  # Initialize with empty tree
  tree.arcs <- matrix(ncol = 3)[-1, ]
  
  stages <- 0  # initialize counter
  stages.arcs <- c()  # vector to store stage number in wich each arc was added
  
  # Iterate while there are more than one component
  while(length(unique(components[, 2])) > 1) {
    
    k <- 1
    nComp <- components[, 2]
    
    # Check each component
    while(k < length(nComp)) {
      
      # Set of component arcs
      sArcs <- matrix(ncol = 3)[-1, ]
      # Set of component nodes
      sNodes <- components[components[, 2] %in% components[k, 2], 1]
      # Review nodes for each component
      for (i in sNodes) {
        
        # Get minimum cost arc to a node from outside the component
        jNodes <- components[-which(components[, 1] %in% sNodes), 1]
        validArcs <- matrix(arcs[arcs[, 1] %in% i & arcs[, 2] 
                                 %in% jNodes, ], ncol = 3)
        if (length(validArcs) > 0) {
          min.arc <- matrix(validArcs[validArcs[, 3] == min(validArcs[, 3]), ],
                          ncol = 3)[1, ]
          # Add arc to S
          sArcs <- rbind(sArcs, min.arc)
        }
        
      }
      
      # Add to T minimum cost arc from S
      ijArc <- matrix(sArcs[sArcs[, 3] == min(sArcs[, 3]), ], ncol = 3)[1, ]
      tree.arcs <- rbind(tree.arcs, ijArc)
      # Check components of the two nodes of selected arc
      iComp <- components[components[, 1] == ijArc[1], 2]
      jComp <- components[components[, 1] == ijArc[2], 2]
      # Merge components
      components[components[, 2] == jComp, 2] <- iComp
      components[components[, 2] > jComp, 2] <- 
        components[components[, 2] > jComp, 2] - 1
      
      # Check next component
      k <- k + 1
      nComp <- components[, 2]
      
    }
    
    # Save in which stage an arc was added to the tree and update
    stages <- stages + 1  # counter
    stages.arcs <- c(stages.arcs,
                     rep(stages, nrow(tree.arcs) - length(stages.arcs)))
    
  }
  
  # Remove component column
  tree.arcs <- matrix(tree.arcs[, -4], ncol = 3)
  # Column names
  colnames(tree.arcs) <- c("ept1", "ept2", "weight")
  
  output <- list("tree.nodes" = nodes, "tree.arcs" = tree.arcs,
                 "stages" = stages, "stages.arcs" = stages.arcs)
  return(output)
  
}
#-----------------------------------------------------------------------------#