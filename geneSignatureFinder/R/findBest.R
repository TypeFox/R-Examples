findBest <-
function(cutoff, ssignature, ccandidate) {
  repeat {
    K <-  length(ccandidate)
    leaveOneOut <- rep(0, K)
    for(k in 1:K) 
      leaveOneOut[k] <- survdiff(stData ~ classify(geData[, c(ssignature, ccandidate[-k])])$clusters)$chisq

    mmax <- max(leaveOneOut)
    if(sum(leaveOneOut == mmax) == K) 
      return(ccandidate)
    
    ccandidate <- ccandidate[leaveOneOut < mmax]
    tmp <-  survdiff(stData ~ classify(geData[, c(ssignature, ccandidate)])$clusters)$chisq
    if(length(ccandidate) < 2)  
      if(tmp > cutoff)
        return(ccandidate) else return(NULL)
  }
}
