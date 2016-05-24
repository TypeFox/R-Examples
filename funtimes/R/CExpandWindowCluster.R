CExpandWindowCluster <- function(e, Euncl){
  if(is.null(dim(Euncl)[2])) {Euncl <- matrix(Euncl, ncol=1)}
  ClusterInclude <- rep(FALSE, dim(Euncl)[2])
  if(length(ClusterInclude)==0) {return(ClusterInclude)}
  if(any(e)) { #Perform it only if we have neighbors
    ClusterInclude[e] <- TRUE
    NeighID <- which(e)
    for (i in 1:sum(e)) {
      if (all(ClusterInclude)) {return(ClusterInclude)}
      tmp <- CExpandWindowCluster(Euncl[!ClusterInclude,NeighID[i]], Euncl[!ClusterInclude,!ClusterInclude])
      ClusterInclude[!ClusterInclude][tmp] <- TRUE
    }
  }  
  return(ClusterInclude)
}