########################################################################################
##computes graph-theoretic information about shortest paths from given source vertices##
##to all target vertices using an adjacency matrix; NOTE: this is an INTERNAL FUNCTION##
##not exported by the package, and as such, it does not provide checks of its argument##
########################################################################################
shortest.paths.information <- function(M){ 
  cc <- .Call("shortestPathsInformation", M)
  
  weight.distances <- cc[[1]] # matrix of the weight-based lengths of the shortest paths from source vertices
  # (row stimuli) to target vertices (column stimuli) (in Fechnerian scaling context,
  # matrices of the oriented Fechnerian distances of the first and second kind)
  predecessors <- cc[[2]]+1  # matrix of the predecessors of the column stimuli in shortest paths from the row stimuli
  # (as source vertices) to the column stimuli (as target vertices)
  # predecessors is an index matrix; in C indices start with 0, thus '+1'
  edge.distances <- matrix(nrow = dim(M)[1], ncol = dim(M)[1]) # matrix of the edge/link based (graph-theoretic) lengths of the shortest paths
  # from source vertices (row stimuli) to target vertices (column stimuli)
  
  for(id in 1:dim(M)[1]){
    lvl <- rep(NA, dim(M)[1])
    lvl[id] <- 0
    distance.2 <- 0
    done.2 <- NA
    used <- numeric()
    while(any(is.na(lvl))){
      done.2 <- which(!is.na(lvl))[!is.element(which(!is.na(lvl)), used)]
      distance.2 <- (distance.2 + 1)
      lvl[which(is.element(predecessors[id,], done.2))] <- distance.2
      used <- append(used, done.2[!is.element(done.2, used)])
    }
    edge.distances[id, ] <- lvl
  }
  
  dimnames(weight.distances) <- dimnames(edge.distances) <- dimnames(predecessors) <- dimnames(M)
  return(list(weight.distances = weight.distances, edge.distances = edge.distances, predecessors = predecessors))
}
