#####################################################
#
#  This file contains the following functions
#
#  FindCG(): find the connected subgraphs with given length
#
#  FindAllCG(): find the connected subgraphs no longer than a given length
#
####################################################

FindCG <- function(adjacency.matrix,cg.initial) {
  ###############################################################################
  # Find all the connected subgraphs with a certain number of nodes 
  # in an undirected graph.
  #
  # Args:
  #   adjacency.matrix: the adjacency matrix of an undirectd graph
  #   cg.initial: the list of all the connected subgraphs with ii nodes
  #
  # Returns: 
  #   cg.new: a matrix with ii+1 colomns, each row of which is 
  #          a connected subgraphs with ii+1 nodes
  ###############################################################################

  if(isSymmetric(adjacency.matrix) == 0){
    stop("The adjacency matrix is not symmetric!")
  }
  p<-dim(adjacency.matrix)[1]
  #set the diagonal of the adjacency matrix to be 0
  if (sum(abs(diag(adjacency.matrix))) > 0) {
    adjacency.matrix <- adjacency.matrix - diag(diag(adjacency.matrix))
  }
  cg.new <- NULL
  #if ((length(dim(cg.initial)) == 0) || (min(dim(cg.initial) == 1))) {
  if (identical(sort(cg.initial),c(1:p)) == 1) {
    cg.initial <- matrix(cg.initial,ncol=1)
  }
  n.cg.ii <- dim(cg.initial)
  if (length(n.cg.ii) == 0) {
    # If there is only 1 cg with ii nodes in a graph, 
    # there is no cg with ii+1 nodes for sure
    stop(paste(paste('No connected subgraphs with ', length(cg.initial)+1), ' nodes detected.'))
  } else{
    n.cg <- n.cg.ii[1]
    ii <- n.cg.ii[2]
    nk <- ii + 1
    for (j in 1:n.cg) {
      # want to construct the connected subgraphs with ii+1 nodes by adding
      # one node to those with ii nodes 
      cg0 <- cg.initial[j, ]
      # list the neighbors of the nodes in cg0
      neighbor.candidate <- which(matrix(adjacency.matrix[cg0, ],nrow = ii) != 0, arr.ind = T)
      neighbor <- neighbor.candidate[, 2]
      neighbor <- unique(setdiff(neighbor, cg0))
      length.cg0 <- length(neighbor)
      if (length.cg0 > 0){
        # create a matrix, each row is a connected subgraph with ii+1 nodes
        # created by adding one node to a connected subgraph with ii nodes
        cg.new.j <- matrix(c(rep(cg0, each=length.cg0), neighbor), nrow=length.cg0)
        cg.new <- rbind(cg.new, cg.new.j)
      }
    }
    if (length(cg.new) == 0) {
      stop(paste(paste('No connected subgraphs with ', nk), ' nodes detected.'))
    } else {
      # sorted the indices in each connected subgraph in ascending order
      cg.new <- t(apply(cg.new, 1, sort))
      # delete the repetitions in the list
      iii <- 1
      nr.cg <- dim(cg.new)[1]
      nl.cg <- dim(cg.new)[2]   
      while (iii < nr.cg){
        # for each current conneced subgraph, 
        # delete the identical connected subgraphs in the later rows. 
        cg.current <- matrix(rep(cg.new[iii, ], (nr.cg - iii)), nrow=nr.cg-iii, byrow=TRUE)
        ix.identical <- rowSums(cg.new[(iii + 1):nr.cg, ] == cg.current)
        id.current <- which(ix.identical == nl.cg) + iii
        diff.current <- setdiff(1:nr.cg, id.current)
        cg.new <- cg.new[diff.current, ]
        nr.cg <- length(diff.current)
        iii <- iii + 1
      }
    }
  }
  return(cg.new)
}

FindAllCG <- function(adjacency.matrix, lc) {
  ###############################################################################
  # Use the function FindCG iteratively to find all the connected subgraphs
  # with no more than lc nodes. 
  #
  # Args:
  #   adjacency.matrix: a graph with p nodes, in the form of adjacency matrix
  # lc: the maximum size of the connected subgraphs to be listed.
  #
  # Returns:
  #   cg.all: a length lc list, whose ii'th element is a matrix whose rows store
  #         the connected subgraphs with ii nodes. 
  ###############################################################################
  if(isSymmetric(adjacency.matrix) == 0){
    stop("The adjacency matrix is not symmetric!")
  }
  p <- dim(adjacency.matrix)[1]
  cg.all <- vector("list", lc)
  # the connected subgraphs with one node are the nodes themselves
  cg.all[[1]] <- 1:p
  if (lc >= 2) {
    support.adjacency <- which(adjacency.matrix != 0, arr.ind = T);
    Ix.support <- (support.adjacency[, 1] < support.adjacency[, 2]);
    cg.all[[2]] <- support.adjacency[Ix.support, ];
  }
  if (lc >= 3) {
    for (ii in 3:lc) {
      ###find all the CG with ii nodes in the graph
      cg.all[[ii]] <- FindCG(adjacency.matrix, cg.all[[ii - 1]]);
    }
  }   
  return(cg.all)
}


