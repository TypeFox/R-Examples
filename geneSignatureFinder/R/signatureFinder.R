signatureFinder <-
function(seedGene, 
  logFilePrefix = "",
  coeffMissingAllowed = 0.75,
  subsetToUse = 1:ncol(geData),
  cpuCluster = NULL, stopCpuCluster = TRUE) {
     
  if(areDataNotLoaded()) return(NULL)
  
  zero = NULL
  
  if(is.null(cpuCluster)) {
    ans <- sequentialSignatureFinder(seedGene, subsetToUse = subsetToUse,
      logFileName = logFilePrefix, coeffMissingAllowed = coeffMissingAllowed)
  } else {
  # export the data to the clusetr of cpu's
  clusterExport(cpuCluster, "geData")
  clusterExport(cpuCluster, "stData")
  # export the funtions to the clusetr of cpu's
  clusterEvalQ(cpuCluster,  library(geneSignatureFinder))
  #instantiate the parallel version of the searching algorithm.
  ans <- parallelSignatureFinder(cpuCluster, seedGene,
      logFileName = logFilePrefix, coeffMissingAllowed = coeffMissingAllowed,
      subsetToUse = subsetToUse, zero = zero)
    if(stopCpuCluster) stopCluster(cpuCluster)
  }
   
  class(ans) <- "gSignature"
  return(ans)
}
