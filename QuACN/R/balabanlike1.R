balabanlike1 <- function(g, dist=NULL){
  # check if g is a graphNEL object
  if(class(g)[1]!="graphNEL"){
    stop("'g' must be a 'graphNEL' object")
  }
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  if(numNodes(g)<3){
    warning("Graps with |V| < 3 result in: Inf! -> value was set to -999")
    return (-999)
  }
  
  mue <- rowSums(dist)
  sigma <- apply(dist,1,max)
  nam <- nodes(g)
  Sj <- lapply(nam,function(n){
    table(dist[n,],exclude=0)
  })
  names(Sj) <- nam
  u <- sapply(nam,function(vi){
    tmp <- sapply(1:sigma[vi],function(j){
      (j*Sj[[vi]][j])/(mue[vi])*log2(j/mue[vi])
    })
    return ((-1)*sum(tmp))
  })

  edges <- .edgesNoDupls(g)
  Ji<-sapply(1:length(edges), function(i){
    start <- names(edges[i])
    targets <- names(edges[[i]])
    1/sqrt(u[start]*u[targets])
  })

  nE <- numEdges(g)
  nV <- numNodes(g)

  return ((nE/(nE-nV+2)) * sum(unlist(Ji)))
}
