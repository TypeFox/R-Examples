"var.xyz" <- function(xyz, weights=TRUE) {
  ## Calculate pairwise distances
  natoms <- ncol(xyz) / 3
  all <- array(0, dim=c(natoms,natoms,nrow(xyz)))
  for( i in 1:nrow(xyz) ) {
    dists <- dist.xyz(xyz[i,])
    all[,,i] <- dists
  }
  
  ## Calculate variance of pairwise distances
  all.vars <- apply(all, 1:2, var)

  if(weights) {
    ## Make the final weights
    wts <- 1 - (all.vars / max(all.vars, na.rm=TRUE))
    wts[is.na(wts)] <- 1
    return(wts)
  }
  else {
    return(all.vars)
  }
}

"var.pdbs" <- function(pdbs, ...) {
  xyz <- pdbs$xyz
  return(var.xyz(xyz, ...))
}

