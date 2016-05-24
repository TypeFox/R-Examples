mahal <- function(q, vc){
  storage <- matrix(NA, nrow(q), nrow(q))
  q <- as.matrix(q)
  ## generate a matrix of Mahalanobis distances between all possible pairs
  for(i in 1:nrow(q)){
    storage[row(q), i] <- sqrt(mahalanobis(x = q, center = q[i,], cov = vc))
  }
  return(storage)
}
