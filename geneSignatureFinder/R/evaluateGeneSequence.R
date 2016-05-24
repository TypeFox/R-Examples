evaluateGeneSequence <-
function(sequenceIDs, 
           coeffMissingAllowed = 0.75) {
    
    n <- nrow(geData)
    m <- ncol(geData)

    result <- list()
    result$startingSignature <-  colnames(geData)[sequenceIDs]
    result$coeffMissingAllowed <- coeffMissingAllowed
    
    if(length(sequenceIDs) > 1) 
      notMissing <- apply(!is.na(geData[, sequenceIDs]), 1, sum) else
    notMissing <- !is.na(geData[, sequenceIDs]) + 0
    notMissing <- notMissing > 0
    
    clusters <- rep(NA, n)
    clusters[notMissing] <- classify(geData[notMissing, sequenceIDs])$clusters
    
    tmp1 <- min(survfit(stData[clusters[notMissing] == 1]~ 1)$surv)
    tmp2 <- min(survfit(stData[clusters[notMissing] == 2]~ 1)$surv)
    if(tmp1 > tmp2) {  
      clusters[notMissing][clusters[notMissing] == 1] <- 0
      clusters[notMissing][clusters[notMissing] == 2] <- 1
    } else clusters[notMissing][clusters[notMissing] == 2] <- 0
    
    result$startingTValue <-  survdiff(stData[notMissing] ~ clusters)$chisq
    result$startingPValue <- 1 - pchisq(result$startingTValue, df = 1)
    
    result$signature <-  result$startingSignature
    result$tValue <-   result$startingTValue 
    result$pValue <-   result$startingPValue
    result$signatureIDs <- sequenceIDs
    names(result$signatureIDs) <- result$signature
    result$classification <- as.factor(clusters)
    levels(result$classification) <- c("good", "poor")
    return(result)
  }
