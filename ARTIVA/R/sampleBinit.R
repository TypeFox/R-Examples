sampleBinit <-
function(Si, sig2, delta2, X, q){
  ### INPUT: s=S[i,],sig2=Sig2[i],delta2,
  ###        X the observed data for predictors.
  ###        q number of predictors
  ### OUTPUT: vector newB.
  ### depends on: q the number of predictors.
  newB <- array(0,q+1)
  
  for(l in which(Si == 1)){
     newB[l] <- rnorm(1, mean=0, sd=sqrt(delta2 * sig2 * t(X[,l]) %*% X[,l]))
  }
  return(newB)
}
