parallelFindBest <-
function(cpuCluster, cutoff, ssignature, ccandidate) {
  repeat {
    K <- length(ccandidate)
    ndx <- matrix(0, ncol = length(ssignature) + K - 1, nrow = K) 
    for(k in 1:K)
      ndx[k,] <- c(ssignature, ccandidate[-k])
    
    clusters <- t(parApply(cpuCluster, ndx, 1, function(nndx)
                           classify(geData[, nndx])$clusters))
    leaveOneOut <- parApply(cpuCluster, clusters, 1, function(cluster)
                            survdiff(stData ~ cluster)$chisq)

    mmax <- max(leaveOneOut)
    if(sum(leaveOneOut == mmax) == K) 
      return(ccandidate)
    
    ccandidate <- ccandidate[leaveOneOut < mmax]
    tmp <-  survdiff(stData ~ classify(geData[, c(ssignature, ccandidate)])$clusters)$chisq
#    if(tmp > cutoff & length(ccandidate) < 2)  
#        return(ccandidate)    
 if(length(ccandidate) < 2)  
      if(tmp > cutoff)
        return(ccandidate) else return(NULL)
  }
}
