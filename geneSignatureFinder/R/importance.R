importance <-
function(aSignatureFinder, deep = FALSE, cpuCluster = NULL, stopCpuCluster = TRUE) {

 if(areDataNotLoaded()) return(NULL)

  if(length(aSignatureFinder$signatureIDs) == 1) {
    message(paste("The importances on the signature starting from", aSignatureFinder$startingSignature, "cannot be computed."))
    return(aSignatureFinder)
  }
  
  signature <- aSignatureFinder$signature
  signatureIDs <- aSignatureFinder$signatureIDs
  tValue <- aSignatureFinder$tValue
 

  ######################
  K <- length(signature)
  
  notMissing <- apply(!is.na(geData[, signatureIDs]), 1, sum)
  notMissing <- notMissing > floor((K-1)^aSignatureFinder$coeffMissingAllowed)

  if(is.null(cpuCluster)) {
    leaveOneOut <- rep(NA, K)
    for(k in 1:K) {
      cluster <- classify(geData[notMissing, signature[-k]])$clusters
      leaveOneOut[k] <- survdiff(stData[notMissing] ~ cluster)$chisq
    }
  } else {
    # export the data to the clusetr of cpu's
    clusterExport(cpuCluster, "geData")
    clusterExport(cpuCluster, "stData")
    # export the funtions to the clusetr of cpu's
    clusterEvalQ(cpuCluster,  library(geneSignatureFinder))
    
    ndx <- matrix(0, ncol = K - 1, nrow = K)    
    for(k in 1:K) ndx[k,] <- signature[-k]
    
    clusters <- t(parApply(cpuCluster, ndx, 1, function(nndx, nnotMissing)
                           classify(geData[nnotMissing, nndx])$clusters, notMissing))
    leaveOneOut <- parApply(cpuCluster, clusters, 1, function(cluster, nnotMissing)
                            survdiff(stData[nnotMissing] ~ cluster)$chisq, notMissing)
  }
  names(leaveOneOut) <-  signature
  aSignatureFinder$L1GeneOutTV <- leaveOneOut
  importance <-  1 - leaveOneOut/tValue
  aSignatureFinder$importance <- importance
#  if(!deep)
  
  if(stopCpuCluster && !is.null(cpuCluster)) stopCluster(cpuCluster)
  return(aSignatureFinder)
  

  
#  oorder <- order(importance)
#  importance <- importance[oorder]
#  signatureIDs <- signatureIDs[oorder]
    
#  res <- rep(0, K)
#  names(res) <- names(signatureIDs)
#  RUNS <- 10*K

#  run <- 1
#  while(run <= RUNS) {
#    print(randomPermutation <- sample(1:K))
#    localSignatureIDs <- signatureIDs[randomPermutation]
#    ans <- rep(0, K)
#    names(ans) <- names(localSignatureIDs)         
#    ans[1] <- importance[names(ans)[1]]
#    ans[1] <- ifelse(ans[1] < 1.1103e-16, 0.0, ans[1])
#    for(k in 2:(K-1)) {
      #print(c(run,k))
#      toDelete <- 1:k
      # qua ci dovrebbe stare il ricalcolo dei missing
#      clustering <- classify(geData[notMissing, localSignatureIDs[-toDelete]])$clusters
#      ans[k] <- 1 - survdiff(stData[notMissing] ~ clustering)$chisq/tValue
#      ans[k] <- ifelse(ans[k] < 1.1103e-16, 0.0, ans[k])
#      ans[k] <- log(ans[k]/aSignatureFinder$actualPValue)/log(10)
#      ans[k] <- tmp <- ifelse(ans[k] <= 0.0, 0.0, ans[k])
#      weight <- ans[1:k]/sum(ans[1:k])
#      weight[which(is.nan(weight))] <- 0
#      ans[1:k] <-  c(ans[1:(k-1)] * (k-1):1, 0.0)
#      ans[1:k] <- (ans[1:k] + weight * tmp) / k:1    
#    }
#    ans[which(ans < 0)] <- 0
#    print(ans)
#    res[randomPermutation] <- res[randomPermutation] + ans
#    run <- run + 1
#  }
#  aSignatureFinder$orderInvariantGeneImportance <- res[aSignatureFinder$signature]
#  return(aSignatureFinder)
}
