## function to compute stress values based on random dissimilarities

randomstress <- function(n, ndim, nrep = 100, type = c("ratio", "interval", "ordinal","mspline")) {

## n ... number of objects
## ndim ... number of dimensions
## nrep ... number of replications  
## type ... MDS type  
  type <- match.arg(type, c("ratio", "interval", "ordinal","mspline"), several.ok = FALSE)
  nn <- n*(n-1)/2                       ## total number of dissimilarities
  
  delta1 <- matrix(NA, n, n)            ## initialize objects
  stressvec <- NULL
  
  for (i in 1:nrep){
    vals <- runif(nn)
    delta1[lower.tri(delta1)] <- vals
    delta <- as.dist(delta1)
    stressvec[i] <- mds(delta, ndim = ndim, type = type)$stress
  }
  return(stressvec)
}