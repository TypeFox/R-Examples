## Hamming distance
##----------------------------------------
hamming <- function(object,...){
  
  adjlist <- object$G
  
  ## for weighted networks, weights must be in [0,1]
  if (object$tag == "undir"){
    dist <- ham.undir(adjlist, object$N, ...)
  } else{
    dist <- ham.dir(adjlist, object$N, ...)
  }
  return(dist)
}


## Useful function for computing Hamming distance
## --------------------------------------------------

## Hamming distance for undirected graph
ham.undir <- function(adjlist, n, ...){
  if (length(adjlist) == 2){
    dist <- sum(abs(adjlist[[1]]-adjlist[[2]]))/(n*(n-1))
    names(dist) <- "H"
  } else {
    idx <- combn(length(adjlist),2)
    tmpdist <- sapply(1:dim(idx)[2],function(x, adjlist, idx){
      sum(abs(adjlist[[idx[1,x]]]-adjlist[[idx[2,x]]]))/(n*(n-1))
    }, adjlist=adjlist, idx=idx)
    dist <- matrix(NA,ncol=length(adjlist), nrow=length(adjlist))
    dist[t(idx)] <- dist[t(idx)[,c(2,1)]] <- tmpdist
    diag(dist) <- 0
  }
  return(dist)
}

## Hamming distance for directed graph
ham.dir <- function(adjlist, n, ...){
  if (length(adjlist) == 2){
    dist <- sum(abs(adjlist[[1]]-adjlist[[2]]))/(2*n*(n-1))
    names(dist) <- "H"
  } else {
    idx <- combn(length(adjlist),2)
    tmpdist <- sapply(1:dim(idx)[2],function(x, adjlist, idx){
      sum(abs(adjlist[[idx[1,x]]]-adjlist[[idx[2,x]]]))/(2*n*(n-1))
    }, adjlist=adjlist, idx=idx)
    dist <- matrix(NA,ncol=length(adjlist), nrow=length(adjlist))
    dist[t(idx)] <- dist[t(idx)[,c(2,1)]] <- tmpdist
    diag(dist) <- 0
  }
  return(dist)
}
