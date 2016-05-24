yhat_SVM <- function(dfTrain, dfTest) {
  Xtr <- as.matrix.data.frame(dfTrain)
  Xte <- as.matrix.data.frame(dfTest)
  p <- ncol(Xtr)-1
  ytr <- Xtr[, p+1]
  Xtr <- Xtr[,-(p+1)]
  Xte <- Xte[,-(p+1)]
  ans <- svm(x=Xtr, y=ytr)
  predict(ans, newdata=Xte)
}
