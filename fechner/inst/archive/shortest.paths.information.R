########################################################################################
##computes graph-theoretic information about shortest paths from given source vertices##
##to all target vertices using an adjacency matrix; NOTE: this is an INTERNAL FUNCTION##
##not exported by the package, and as such, it does not provide checks of its argument##
########################################################################################
shortest.paths.information <-
  function(M){
  
    # M: an adjacency matrix (in Fechnerian scaling context, matrices of the psychometric increments
    #    of the first and second kind)
  
    n <- dim(M)[1]
    weight.distances <- matrix(nrow = n, ncol = n)  # matrix of the weight-based lengths of the shortest paths from source vertices
    # (row stimuli) to target vertices (column stimuli) (in Fechnerian scaling context,
    # matrices of the oriented Fechnerian distances of the first and second kind)
    edge.distances <- matrix(nrow = n, ncol = n)  # matrix of the edge/link based (graph-theoretic) lengths of the shortest paths
    # from source vertices (row stimuli) to target vertices (column stimuli)
    predecessors <- matrix(nrow = n, ncol = n)  # matrix of the predecessors of the column stimuli in shortest paths from the row stimuli
    # (as source vertices) to the column stimuli (as target vertices)
    
    for(id in 1:n){
      node.from <- id  # node.from: a given source vertex (row stimulus) for which to determine information about shortest paths
      # to the column stimuli (as target vertices)
      distance <- rep(Inf, n)
      lvl <- rep(NA, n)
      pred <- rep(0, n)
      done <- rep(FALSE, n)
      distance[node.from] <- 0
      lvl[node.from] <- 0
      for(i in 1:n){
        node.closest <- (-1)
        min.dist <- Inf
        for(j in 1:n){
          if(!done[j]){
            if(distance[j] <= min.dist){
              min.dist <- distance[j]
              node.closest <- j
            }
          }
        }
        done[node.closest] <- TRUE
        for(j in 1:n){
          if(!done[j]){
            if((distance[node.closest] + M[node.closest, j]) < distance[j]){
              distance[j] <- (distance[node.closest] + M[node.closest, j])
              pred[j] <- node.closest
            }
          }
        }
      }
      distance.2 <- 0
      done.2 <- NA
      used <- numeric()
      while(any(is.na(lvl))){
        done.2 <- which(!is.na(lvl))[!is.element(which(!is.na(lvl)), used)]
        distance.2 <- (distance.2 + 1)
        lvl[which(is.element(pred, done.2))] <- distance.2
        used <- append(used, done.2[!is.element(done.2, used)])
      }
      weight.distances[id, ] <- distance
      edge.distances[id, ] <- lvl
      predecessors[id, ] <- pred
    }
    dimnames(weight.distances) <- dimnames(edge.distances) <- dimnames(predecessors) <- dimnames(M)
    
    return(list(weight.distances = weight.distances, edge.distances = edge.distances, predecessors = predecessors))
  }
