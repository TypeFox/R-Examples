# ************************************************************************
# Predict the class of new data using an existing ensemble
# ************************************************************************

cfPredict <- function(ensObj, newInputData) {
  predMatr = finalClasses = confScores = c()
  
  # If missing new data or null, throw an error and stop
  if (missing(newInputData) || is.null(newInputData)) {
    stop("Please provide a matrix with test data.")
  } else {
    
    # Convert the input arguments into the right format
    newData <- as.matrix(as.data.frame(newInputData))
    ensNum  <- length(ensObj$svmModel)
    
    # Loop through each SVM model in the ensemble and predict the new data
    for (i in 1:ensNum) {
      newModel <- ensObj$svmModel[[i]]
      predTest <- predict(newModel, newData)
      predMatr <- cbind(predMatr, as.character(predTest))
    }
    
    colnames(predMatr) <- paste("model", 1:ensNum, sep="")
    
    if (is.null(rownames(predMatr))) {
      rownames(predMatr) <- paste("sample", 1:nrow(predMatr), sep="")
    }
    
    confList <- lapply(apply(predMatr, 1, function(x) list(table(x))), "[[", 1)
    
    for (j in 1:length(confList)) {
      finalClasses <- c(finalClasses, names(which.max(confList[[j]])))
      confScores   <- c(confScores, max(prop.table(confList[[j]]) * 100))
    }
    
    # Combine the results into a data frame
    finalMatr <- data.frame(finalClasses, as.numeric(confScores))
    rownames(finalMatr) <- rownames(newData)
    colnames(finalMatr) <- c("Voted Class", "Conf Score(%)")
    
    return (finalMatr)
  }
}
