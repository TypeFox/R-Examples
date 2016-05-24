Trunkpoint <-
function(B, K, pvalues) {
  trunks <- matrix(0,ncol=K, nrow=B+1)
  trunks[,1] <- -log(pvalues[,1])
  if(K==1) return(trunks)
  for(j in 2:K) {
    trunks[,j] <- trunks[,j-1]+(-log(pvalues[,j]))
  }			
  return(trunks)
}
