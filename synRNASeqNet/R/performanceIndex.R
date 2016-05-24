performanceIndex <-
function(testNet, gsNet){
  ans <- as.matrix(sort(unique(as.vector(testNet)), decreasing = TRUE))
  ans <- cbind(ans, matrix(as.numeric(sapply(ans[, 1], function(x) performanceNET(testNet >= x, gsNet))), ncol = 11, byrow = T))
  colnames(ans) <- c("Thresh", "PredPairs", "ValidPairs", "TP", "FP", "FN", "TN", "Recall", "FPR", "Precision", "Accuracy", "Fscore")
  ans <- as.data.frame(ans)
  
  return(ans)
}

