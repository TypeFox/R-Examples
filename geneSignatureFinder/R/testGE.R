testGE <-
function(aSignatureFinder, permutationReplications  = 1000, cpuCluster = NULL, stopCpuCluster = TRUE) {
  
   if(areDataNotLoaded()) return(NULL)
  
  signatureIDs <- aSignatureFinder$signatureIDs
  L <- length(signatureIDs)
  n <- nrow(geData)
  
  if(!is.null(cpuCluster)) {
    clusterExport(cpuCluster, "geData")
    permutations <- matrix(0, nrow = n, ncol = permutationReplications)
    for(b in 1:permutationReplications)
      permutations[, b] <- sample(aSignatureFinder$classification)
    tmp <- parApply(cpuCluster, permutations, 2, 
                    function(randomSample, LL, ids) {
                      ans <- rep(0, LL)
                      for(l in 1:LL)
                        ans[l] <- t.test(geData[randomSample == 1, ids[l]], 
                                         geData[randomSample == 2, ids[l]], 
                                         var.equal = FALSE)$statistic
                      return(ans)
                    }, L, signatureIDs)
    if(stopCpuCluster) stopCluster(cpuCluster)
  } else {
    tmp <- matrix(0, nrow = L, ncol = permutationReplications)
    for(b in 1:permutationReplications) {
      randomSample <- sample(unclass(aSignatureFinder$classification))
      for(l in 1:L)
        tmp[l, b] <- t.test(geData[randomSample == 1, signatureIDs[l]], 
                            geData[randomSample == 2, signatureIDs[l]], 
                            var.equal = FALSE, na.rm = TRUE)$statistic
    }
  }
  tValue <- rep(0, L)
  meanDifference <- rep(0, L)
  medianDifference <- rep(0, L)
  groupMean <- matrix(0, ncol = 2, nrow = L)
  groupMedian <- matrix(0, ncol = 2, nrow = L)
  for(l in 1:L) {
    tValue[l] <- t.test(geData[aSignatureFinder$classification == "good", signatureIDs[l]], 
                        geData[aSignatureFinder$classification == "poor", signatureIDs[l]], 
                        var.equal = FALSE, na.rm = TRUE)$statistic
    groupMean[l,] <- c(mean(geData[aSignatureFinder$classification == "good", signatureIDs[l]], na.rm = TRUE),
                       mean(geData[aSignatureFinder$classification == "poor", signatureIDs[l]], na.rm = TRUE))
    
    groupMedian[l,] <- c(median(geData[aSignatureFinder$classification == "good", signatureIDs[l]], na.rm = TRUE),
                       median(geData[aSignatureFinder$classification == "poor", signatureIDs[l]], na.rm = TRUE))
    meanDifference[l] <- groupMean[l,1] - groupMean[l,2]
    medianDifference[l] <- groupMedian[l,1] - groupMedian[l,2]
    }
  rownames(groupMean) <- aSignatureFinder$signature
  colnames(groupMean) <- c("good", "poor")
  rownames(groupMedian) <- aSignatureFinder$signature
  colnames(groupMedian) <- c("good", "poor")
  pValue <- rep(0, L)
  for(l in 1:L)
    pValue[l] <- ifelse(tValue[l] > 0, mean(tValue[l] < tmp[l,]), mean(tValue[l] > tmp[l,]))
  
  aSignatureFinder$groupMedian <- groupMedian
  aSignatureFinder$medianAbsDifference <- abs(medianDifference)
  names(aSignatureFinder$medianAbsDifference) <- aSignatureFinder$signature
  aSignatureFinder$groupMean <- groupMean
  aSignatureFinder$meanAbsDifference <- abs(meanDifference)
  names(aSignatureFinder$meanAbsDifference) <- aSignatureFinder$signature
  aSignatureFinder$meanDifferenceTValue <- tValue
  names(aSignatureFinder$meanDifferenceTValue) <- aSignatureFinder$signature
  aSignatureFinder$meanDifferencePValue <- pValue
  names(aSignatureFinder$meanDifferencePValue) <- aSignatureFinder$signature
  return(aSignatureFinder)
}
