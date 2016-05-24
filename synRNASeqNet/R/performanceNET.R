performanceNET <-
function(predAdjMat, valAdjMat){
  TP <- sum(predAdjMat != 0 & valAdjMat != 0)
  FP <- sum(predAdjMat != 0 & valAdjMat == 0)
  FN <- sum(predAdjMat == 0 & valAdjMat != 0)
  TN <- sum(predAdjMat == 0 & valAdjMat == 0)
  
  nEdges <- TP + FP
  recall <- TP/(TP + FN)
  FPR <- FP/(FP + TN)
  precision <- TP/(TP + FP)
  precision[is.na(precision)] <- 0
  accuracy <- (TP + TN)/(TP + FN + FP + TN)
  Fscore <- 2*(precision*recall)/(precision + recall)
  Fscore[is.na(Fscore)] <- 0
  
  ans <- data.frame(nEdges, sum(valAdjMat != 0), TP, FP, FN, TN, recall, FPR, precision, accuracy, Fscore)
  names(ans) <- c("PredPairs" ,"ValidPairs", "TP", "FP", "FN", "TN", "Recall", "FPR", "Precision", "Accuracy", "Fscore")
  #ans <- round(ans, 4)
  
  return(ans)
}
