topologicalInfoContent <- function(g, dist=NULL, deg=NULL){
  if(class(g)[1]!="graphNEL"){
  stop("'g' must be a 'graphNEL' object")
  }
  stopifnot(.validateGraph(g))
  
  if(is.null(dist)){
    dist <- distanceMatrix(g)
  }
  if(is.null(deg)){
    deg <- graph::degree(g)
  }
  On <- .cardNi(g,dist,deg)
  pis <- On/sum(On)
  In=pis*log2(pis)
  Iorb <- (-1)*sum(In)
  ret <- list()
  ret[["entropy"]] <- Iorb
  ret[["orbits"]] <- On
  return(ret)
}

.getTopology <- function(g,dist,deg){
    nodes <- nodes(g)
    topology<-lapply(nodes,function(vi){
    dist.vi<-dist[vi,]
    max.dist <- max(dist.vi)
    dist.vi <- dist.vi[dist.vi!=0]
    deg.vi <- deg[(names(dist.vi))]
    names(deg.vi) <- dist.vi
    deg.vi <- deg.vi[order(names(deg.vi))]
    tmp.names <-names(deg.vi)
    unames <- unique(names(deg.vi))
    deg.vi <- unlist(lapply(unames,function(un){
      csel <- grep(paste("^",un,"$",sep=""),names(deg.vi))
      return(sort(deg.vi[csel]))
    }))
    return(deg.vi)
  })
  names(topology) <- nodes
  return(topology)
}

.getOrbits <- function(top){
  rest <- names(top)
  i <- 1
  erg <- list()
  while(length(rest)>0){
    rem <- c()
    vi <- rest[1]
    erg.tmp <- vi
    rest <- rest[-1]
    #rest.tmp <- rest
    if(length(rest)>0){
    for(j in length(rest):1){
      if(identical(top[[vi]],top[[rest[j]]])){
        erg.tmp <- cbind(erg.tmp,rest[j])
        rest <- rest[-j]
      }
    }
  }
    erg[[i]]<-sort(erg.tmp)
    i<-i+1
  }
  return(erg)
}

.cardNi <- function(g, dist, deg){
  top <- .getTopology(g,dist,deg)
  orbits <- .getOrbits(top)
  return(sapply(orbits,length))
}
