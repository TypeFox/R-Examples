ensembleTable <-
function(aSearchResults) {
  K <- length(aSearchResults)
  seedGenes <- aSearchResults[[1]]$startingSignature
  for(k in 2:K)
    seedGenes <- c(seedGenes, aSearchResults[[k]]$startingSignature)
  
  ensemble <- aSearchResults[[1]]$signature
  for(k in 2:K) 
    ensemble <- c(ensemble, aSearchResults[[k]]$signature)
  ensemble <- table(ensemble)
  ensemble <- ensemble[rev(order(ensemble))]
  ensembleNames <- names(ensemble)
  nensemble <- length(ensembleNames)

  tableOfImportance <- matrix(0, nrow = nensemble, ncol = K)
  colnames(tableOfImportance) <- seedGenes
  rownames(tableOfImportance) <- ensembleNames
  for(k in 1:K) {
    if(length(aSearchResults[[k]]$signature) == 1) 
      next
    tableOfImportance[aSearchResults[[k]]$signature, k] <- aSearchResults[[k]]$importance
  }
  
  tableOfImportance <- round(10000 * tableOfImportance)/100
  tableOfImportance[which(tableOfImportance == "0")] <- ""
  meanGeneImportance <- rep(0, nensemble)
  names(meanGeneImportance) <- ensembleNames
  for(k in 1:nensemble) {
    gene <- ensembleNames[k]
    numberOfPresences <- 0
    for(kk in 1:K) {
      if(gene %in% names(aSearchResults[[kk]]$importance)) {
        ndx <- which(gene == names(aSearchResults[[kk]]$importance))
        numberOfPresences <- numberOfPresences + length(ndx)      
        meanGeneImportance[k] <- meanGeneImportance[k] + sum(aSearchResults[[kk]]$importance[ndx])
      }
    }
    meanGeneImportance[k] <- meanGeneImportance[k]/numberOfPresences
  }
  correctedMeanGeneImportance <- meanGeneImportance * ensemble/K
#missingPct <- round(100 * apply(geData[, ensembleNames], 2, function(xx) sum(is.na(xx)))/m, 2)
  summaryTable <- cbind(#missingPct, 
    ensemble, 
    round(10000 * meanGeneImportance)/100, 
    round(10000 * correctedMeanGeneImportance)/100)
#colnames(summaryTable) <- c("missing", "counts", "importance", "wImportance")
  colnames(summaryTable) <- c("counts", "importance", "wImportance")
  ensemble <- as.data.frame(cbind(summaryTable, tableOfImportance))
  return(ensemble)
}
