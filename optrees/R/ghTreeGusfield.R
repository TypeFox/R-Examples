#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimum Cut Tree Problems                                                   #
#-----------------------------------------------------------------------------#

# ghTreeGusfield --------------------------------------------------------------
#' Gomory-Hu tree with the Gusfield's algorithm
#'  
#' Given a connected weighted and undirected graph, the \code{ghTreeGusfield}
#' function builds a Gomory-Hu tree with the Gusfield's algorithm.
#' 
#' @details The Gomory-Hu tree was introduced by R. E. Gomory and T. C. Hu in
#' 1961. Given a connected weighted and undirected graph, the Gomory-Hu tree
#' is a weighted tree that contains the minimum s-t cuts for all s-t pairs
#' of nodes in the graph. Gomory and Hu also developed an algorithm to find it
#' that involves maximum flow searchs and nodes contractions.
#' 
#' In 1990, Dan Gusfield proposed a new algorithm that can be used to find a
#' Gomory-Hu tree without nodes contractions and simplifies the implementation.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' 
#' @return \code{ghTreeGusfield} returns a list with:
#' tree.nodes vector containing the nodes of the Gomory-Hu tree.
#' tree.arcs matrix containing the list of arcs of the Gomory-Hu tree.
#' stages number of stages required.
#' 
#' @references R. E. Gomory, T. C. Hu. Multi-terminal network flows. Journal 
#' of the Society for Industrial and Applied Mathematics, vol. 9, 1961.
#' 
#' Dan Gusfield (1990). "Very Simple Methods for All Pairs Network Flow 
#' Analysis". SIAM J. Comput. 19 (1): 143-155.
#' 
#' @seealso A more general function \link{getMinimumCutTree}.

ghTreeGusfield <- function(nodes, arcs) {
  
  # Previous we have a order vector of nodes
  # Start a tree with one node
  nodesT1 <- nodes[1]
  arcsT1 <- matrix(ncol = 4)[-1, ]
  
  # Iterate adding one arc between node i and one node of the tree  
  for (i in 2:length(nodes)) {
    
    # Method to chose one node of the tree
    nodesT <- nodesT1
    arcsT <- arcsT1
    
    # Iterate until have a tree with one node
    while (length(nodesT) > 1) {
      
      # Search a-b arc with minimum weight
      min.arc <- which(arcsT[, 3] == min(arcsT[, 3]))[1]
      # This arc has the weight of the minimum a-b cut in the original graph
      a <- arcsT[min.arc, 1]
      b <- arcsT[min.arc, 2]
      # Remove arc by make it and arc with zero capacity
      arcsT[min.arc, 4] <- 0
      
      # Duplicate and order arcs to find the cut
      arcsT2 <- rbind(arcsT, matrix(c(arcsT[, 2], arcsT[, 1],
                                      arcsT[, 3], arcsT[, 4]), ncol = 4))
      arcsT2 <- arcsT2[order(arcsT2[, 1], arcsT2[, 2]), ]
      # Have two components
      TaTbCut <- findstCut(nodesT, arcsT2, a, b)
      
      # Extract arcs of the two components
      nodesTa <- TaTbCut$s.cut
      arcsTa <- matrix(arcsT[which(arcsT[, 1] %in% nodesTa
                                   & arcsT[, 2] %in% nodesTa), ], ncol = 4)
      nodesTb <- TaTbCut$t.cut
      arcsTb <- matrix(arcsT[which(arcsT[, 1] %in% nodesTb 
                                   & arcsT[, 2] %in% nodesTb), ], ncol = 4)
      # And we have two components in the original graph
      # Use function findMinCut to recover them
      abCut <- findMinCut(nodes, arcs, source.node = a, sink.node = b)
      # Select nodes and arcs connected with node i
      if (i %in% abCut$s.cut) {
        nodesT <- nodesTa
        arcsT <- arcsTa
      } else {
        nodesT <- nodesTb
        arcsT <- arcsTb
      }
      
    }
    
    # At the end we hace one tree with only one node
    nodesT
    # Compute minimum cut i-k
    ikCut <- findMinCut(nodes, arcs, source.node = nodesT, sink.node = i)
    iCut <- ikCut$s.cut
    kCut <- ikCut$t.cut
    ikFlow <- ikCut$max.flow
    # Connect node from tree with i node with weigth equal to minimum cut i-k
    nodesT1 <- c(nodesT1, i)
    arcsT1 <- rbind(arcsT1, c(nodesT, i, ikFlow, ikFlow))
    
  }
  
  # Remove columns of capacities
  tree.arcs <- arcsT1[, -4]
  # Order arcs
  tree.arcs <- tree.arcs[order(tree.arcs[, 1], tree.arcs[, 2]), ]
  # Column names
  colnames(tree.arcs) <- c("ept1", "ept2", "weight")
    
  # Build output
  output <- list("tree.nodes" = nodes, "tree.arcs" = tree.arcs,
                 "stages" = length(nodes))
  return(output)
  
}
#-----------------------------------------------------------------------------#