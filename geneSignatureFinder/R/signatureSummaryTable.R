signatureSummaryTable <-
function(aSignatureFinder) {

    if(length(aSignatureFinder$signature) > 1) {
      n <- nrow(geData)
      missingPct <- round(100 * apply(geData[, aSignatureFinder$signatureIDs], 2, function(xx) sum(is.na(xx)))/n, 2)
      signatureSummary <- cbind(missingPct)  
      rownames(signatureSummary) <- aSignatureFinder$signature
      if(!is.null(aSignatureFinder$importance)) {
        signatureSummary <- cbind(signatureSummary, 
                                  round(aSignatureFinder$importance*100, 2))
        tmp <- length(colnames(signatureSummary))
        colnames(signatureSummary)[tmp] <- "importance"
      }
      if(!is.null(aSignatureFinder$groupMedian)) {
        signatureSummary <- cbind(signatureSummary, 
                                  round(aSignatureFinder$groupMedian, 2))
        tmp <- length(colnames(signatureSummary))
        colnames(signatureSummary)[c(tmp-1, tmp)] <- c("median level in good", "median level in poor")
        
        signatureSummary <- cbind(signatureSummary, 
                                  round(aSignatureFinder$groupMean, 2),
                                  round(aSignatureFinder$meanAbsDifference, 2))
        tmp <- length(colnames(signatureSummary))
        colnames(signatureSummary)[c(tmp-2, tmp-1, tmp)] <- c("mean level in good", 
                                                              "mean level in poor", "mean absolute difference")
        signatureSummary <- cbind(signatureSummary, 
                                  round(aSignatureFinder$meanDifferenceTValue, 3),
                                  aSignatureFinder$meanDifferencePValue)
        tmp <- length(colnames(signatureSummary))
        colnames(signatureSummary)[c(tmp-1, tmp)] <- c("Welch's t-Value", "p-Value")
      }
      signatureSummary <- as.data.frame(signatureSummary)
      return(signatureSummary)
    } else {
      message("The signature has length 1")
      return(NULL)
    }
}
