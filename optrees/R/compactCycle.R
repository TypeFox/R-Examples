#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Arborescence Problems                                               #
#-----------------------------------------------------------------------------#

# compactCycle ----------------------------------------------------------------
#' Compacts the nodes in a cycle into a single node
#'
#' Given a directed graph with a cycle, \code{compactCycle} compacts all the 
#' nodes in the cycle to a single node called super.node. The function uses the
#' first and the last node of the cycle as a fusion point and obtains the costs
#' of the incoming and outgoing arcs of the new node.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param cycle vector with the original nodes in the cycle.
#'
#' @return \code{compactCycle} returns the nodes and the list of arcs forming a
#' new graph with the compressed cycle within a super.node. Also returns a list
#' of the correspondences between the nodes of the new graph and the nodes of 
#' the previous graph.
#' 
#' @seealso This function is an auxiliar function used in 
#' \link{msArborEdmonds} and \link{getMinimumArborescence}.

compactCycle <- function(nodes, arcs, cycle) {
  
  # Previous
  Cmat <- ArcList2Cmat(nodes, arcs)  # work with cost matrix
  
  # Initialize
  super.node <- min(cycle)  # minimun node of the cycle as super.node
  iCosts <- Cmat[super.node, -cycle[2:(length(cycle)-1)]]  # new node file
  jCosts <- Cmat[-cycle[2:(length(cycle)-1)], super.node]  # new node column
  
  # Matches between nodes
  match.nodes <- matrix(nodes, nrow = length(nodes), ncol = 2)
  # First column previous list of nodes, second column new list of nodes
  match.nodes[match(cycle[-length(cycle)], match.nodes[,1]), 2] <- super.node
  # Update nodes with larger label than super.node
  new.nodes <- which(match.nodes[, 2] > super.node)
  match.nodes[new.nodes, 2] <- c((super.node+1):(super.node+length(new.nodes)))
  
  # Fill new file with minimun arcs from removed nodes
  # Costs from leaving arcs
  cost.arcs <- Cmat[cycle[1:(length(cycle)-1)], -cycle[2:(length(cycle)-1)]]
  # Add to fileCost the minimum cost arcs
  for (j in 1:length(iCosts)) {
    # Iterate by columns
    iCosts[j] <- min(cost.arcs[, j])  # save minimum cost arc
  }
  
  # Fill new column with minimum arcs from removed nodes
  # Costs from reaching arcs
  cost.arcs <- Cmat[-cycle[2:(length(cycle)-1)], cycle[1:(length(cycle)-1)]]
  # Add to colCost the minimum cost arcs
  for (i in 1:length(jCosts)) {
    # Iterate by files
    jCosts[i] <- min(cost.arcs[i, ])  # save minimum cost arc
  }
  
  # New vector of nodes
  nodes <- c(1:(length(nodes) - (length(cycle) - 2)))
  # Compact nodes of the cycle in super.node
  Cmat <- Cmat[-cycle[2:(length(cycle)-1)],-cycle[2:(length(cycle)-1)]]
  # Add costs calculated before
  Cmat[super.node, ] <- iCosts
  Cmat[, super.node] <- jCosts
  
  # Return list of arcs
  arcs <- Cmat2ArcList(nodes, Cmat, directed = TRUE)
  
  return(list("nodes"= nodes, "arcs" = arcs, "super.node" = super.node,
              "matches" = match.nodes))
}
#-----------------------------------------------------------------------------#