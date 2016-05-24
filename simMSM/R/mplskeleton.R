mplskeleton <- function(tmat){
  mpl <- vector("list", nrow(tmat))
  names(mpl) <- paste("from", 1:nrow(tmat), sep = ".")
  for(k in 1:nrow(tmat)){
    mpl[[k]] <- vector("list", 4)
    names(mpl[[k]]) <- c("from", "all.to", "bhr", "eta")
    mpl[[k]]$from <- k
    mpl[[k]]$bhr <- mpl[[k]]$eta <- vector("list", ncol(tmat))
    if(any(tmat[k, ])){
      mpl[[k]]$all.to = which(tmat[k, ])
    }
    names(mpl[[k]]$bhr) <- paste("bhr.", k, 1:ncol(tmat), sep = "")
    names(mpl[[k]]$eta) <- paste("eta.", k, 1:ncol(tmat), sep = "")
    for(l in 1:ncol(tmat)){
      if(!tmat[k, l]){
        mpl[[k]]$bhr[[l]] <- function(t){return(0)}
        mpl[[k]]$eta[[l]] <- function(x.i, t){return(0)}
      }
    }
  }
  return(mpl)
}