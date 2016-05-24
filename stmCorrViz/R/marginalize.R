marginalize <-
function(members, beta, weights) {  
  w <- weights[members]/sum(weights[members])
  pvec <- colSums(beta[members,,drop=FALSE]*w)
  words <- order(pvec, decreasing=TRUE)
  return(list(beta=pvec, indices=words))
}
