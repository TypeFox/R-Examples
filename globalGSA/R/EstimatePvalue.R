EstimatePvalue <-
function(B, K, pvalues) {
  WbK <- Trunkpoint(B, K, pvalues)
  sb <- matrix(0,ncol=K, nrow=B+1)
  for(j in 1:K) {
    for(i in 1:(B+1)) {
      sb[i,j] <- sum(WbK[,j] >= WbK[i,j])/(B+1)
    }
  }  
  return(sb)
}
