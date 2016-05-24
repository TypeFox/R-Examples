removeGeneFrom <-
function(aSignatureFinder, cutoff = 0.0) {
  
  if(length(aSignatureFinder$signature) == 1)  {
    message("Signature of length = 1: nothing to do.")
    return(NULL)
  }  

  bby <- "importance"
  if(bby == "importance" && is.null(aSignatureFinder$importance)) {
    message("The importances have to computed before calling this function.")
    return(aSignatureFinder)
  }
  
  toRemove <- which.min(aSignatureFinder$importance)
  if(length(toRemove) > 1)
    toRemove <- sample(toRemove, size = 1)
    
  if(aSignatureFinder$importance[toRemove] <= cutoff) {
    message("Removing ", aSignatureFinder$signature[toRemove], 
            " with importance: ", round(aSignatureFinder$importance[toRemove], 3))
    if(is.null(aSignatureFinder$removedGene)) {
      aSignatureFinder$removedGene <- aSignatureFinder$signature[toRemove]  
      aSignatureFinder$originalSignature <- aSignatureFinder$signature
      aSignatureFinder$originalClassification <- aSignatureFinder$classification
      aSignatureFinder$originalTValue <-  aSignatureFinder$tValue
      aSignatureFinder$originalPValue <- aSignatureFinder$pValue
    } else
        aSignatureFinder$removedGene <- c(aSignatureFinder$removedGene, aSignatureFinder$signature[toRemove])
        
    aSignatureFinder$signatureIDs <- aSignatureFinder$signatureIDs[-toRemove]
    aSignatureFinder$signature <-  aSignatureFinder$signature[-toRemove]
      
    n <- nrow(geData)
    m <- ncol(geData)
    
    if(length(aSignatureFinder$signatureIDs) > 1) 
      notMissing <- apply(!is.na(geData[, aSignatureFinder$signatureIDs]), 1, sum) else
    notMissing <- !is.na(geData[, aSignatureFinder$signatureIDs]) + 0
    notMissing <- notMissing > 0
    
    clusters <- rep(NA, n)
    clusters[notMissing] <- classify(geData[notMissing, aSignatureFinder$signatureIDs])$clusters
 
    
    clusters <- goodAndPoorClassification(clusters) # 08/04/2012
#    tmp1 <- min(survfit(stData[clusters[notMissing] == 1]~ 1)$surv)
#    tmp2 <- min(survfit(stData[clusters[notMissing] == 2]~ 1)$surv)
#    if(tmp1 > tmp2) {  
#      clusters[notMissing][clusters[notMissing] == 1] <- 0
#      clusters[notMissing][clusters[notMissing] == 2] <- 1
#    } else clusters[notMissing][clusters[notMissing] == 2] <- 0
    
#    aSignatureFinder$startingTValue <-  aSignatureFinder$tValue
#    aSignatureFinder$startingPValue <- aSignatureFinder$pValue
    aSignatureFinder$tValue <-  survdiff(stData[notMissing] ~
 clusters[notMissing])$chisq
    aSignatureFinder$pValue <- 1 - pchisq(aSignatureFinder$tValue, df = 1)
    
    aSignatureFinder$classification <- clusters #08/04/2012
    #aSignatureFinder$classification <- as.factor(clusters)
    #levels(aSignatureFinder$classification) <- c("good", "poor")
    if(length(aSignatureFinder$signatureIDs) == 1) {
      message("The signature con no longer be reduced.")
      return(aSignatureFinder)
    }
    aSignatureFinder <- importance(aSignatureFinder)
    
    if(!is.null(aSignatureFinder$groupMedian)) {
      aSignatureFinder$groupMedian <- NULL
      aSignatureFinder$medianAbsDifference <- NULL
      aSignatureFinder$groupMean <- NULL
      aSignatureFinder$meanAbsDifference <- NULL
      aSignatureFinder$meanDifferenceTValue <- NULL      
      aSignatureFinder$meanDifferencePValue <- NULL
    }
    return(aSignatureFinder)
  } else {
    message("No gene to remove at cutoff level: ", cutoff)
    return(NULL)
  }  
}
