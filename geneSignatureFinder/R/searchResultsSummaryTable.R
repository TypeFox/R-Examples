searchResultsSummaryTable <-
function(aSearchResults) {
  K <- length(aSearchResults)
  seedGenes <- aSearchResults[[1]]$signatureName
  for(k in 2:K)
    seedGenes <- c(seedGenes, aSearchResults[[k]]$signatureName)

  tableOfSignatures <- matrix("", ncol = 5, nrow = K)
  rownames(tableOfSignatures) <- seedGenes
  colnames(tableOfSignatures) <- c("length", "tValue", "log(pValue)", "tValue improvement", "signature")
  for (k in 1:K) 
    tableOfSignatures[k, ] <- c(length(aSearchResults[[k]]$signature), 
                                round(aSearchResults[[k]]$tValue, 3), 
                                round(log(aSearchResults[[k]]$pValue)/log(10), 3), 
                                round((1 - aSearchResults[[k]]$startingTValue/aSearchResults[[k]]$tValue) * 100, 2),
                                paste(aSearchResults[[k]]$signature, collapse = ", "))
  tableOfSignatures <- as.data.frame(tableOfSignatures)
  return(tableOfSignatures)
}
