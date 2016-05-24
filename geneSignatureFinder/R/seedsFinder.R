seedsFinder <-
function(cutoff = 1.95,
                        evaluateBICs = TRUE,
                        cpuCluster = NULL) { 
                        
 if(areDataNotLoaded()) return(NULL)
  
  if(!is.null(cpuCluster)) {
    message("Computation started at ", eTime <- Sys.time())
    # export the data to the cluster of cpu's
    clusterExport(cpuCluster, "geData")
    clusterExport(cpuCluster, "stData")
    # export the funtions to the cluster of cpu's
    clusterEvalQ(cpuCluster,  library(geneSignatureFinder))
    #instantiate the parallel version of the searching algorithm.
	  ans <- parApply(cpuCluster, geData, 2, checkGeneForSeed, 
                    cutoff = cutoff, evaluateBICs = evaluateBICs)
    stopCluster(cpuCluster)
    message("Computation ended at ", eTime <- Sys.time())
    return(t(ans))
  }

  n <- nrow(geData)
  m <- ncol(geData)
  ans <- matrix(0, nrow = m, ncol = 4+n)
  colnames(ans) <- c("tValue", "pValue", "bic1", "bic2", rownames(geData))
  rownames(ans) <- colnames(geData)
  j <- 1
  while(j <= m) {
    cat(paste(j,"::", colnames(geData)[j], "\n"))
    ans[j, ] <- checkGeneForSeed(geData[,j],
                                evaluateBICs = evaluateBICs,
                                cutoff = cutoff)
    j <- j + 1
  }
  return(ans)
}
