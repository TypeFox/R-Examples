pmeskeleton <- function(tmat){
  pme <- vector("list", nrow(tmat))
  names(pme) <- paste("from", 1:nrow(tmat), sep = ".")
  for(k in 1:nrow(tmat)){
    pme[[k]] <- vector("list", nrow(tmat))
    names(pme[[k]]) <- paste("eta.", k, 1:nrow(tmat), sep = "")
    for(l in 1:ncol(tmat)){
      if(!tmat[k, l]){
        pme[[k]][[l]] <- function(x){return(0)}
      }
    }
  }
  return(pme)
}