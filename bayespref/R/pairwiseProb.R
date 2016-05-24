pairwiseProb <-
function(prefres=NULL,burn=0){
  pref<-prefres$PopPref
  mcmc<-dim(pref)[2]
  ## remove burnin and collect prefernce data
  pref<-pref[,(burn + 1):mcmc]
  nCat<-dim(pref)[1]
  pairP<-array(dim=c(nCat,nCat))
  for (i in 1:nCat){
    for (j in 1:nCat){
      if (i == j){ ## diagonal of matrix
        pairP[i,j]<-0
      }
      else {
        pairP[i,j]<-sum(pref[i,] > pref[j,]) / length(pref[i,])
      }
    }
  }
  return(pairP)
}

