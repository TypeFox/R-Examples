yhat_lars <-
function(dfTrain, dfTest, normalize=TRUE){
  Xtr <- as.matrix.data.frame(dfTrain)
  ytr <- Xtr[,ncol(Xtr)]
  Xtr <- Xtr[,-ncol(Xtr)]
  ans <- lars(Xtr, ytr, type="lasso", normalize=normalize)
  iopt <- which.min(summary(ans)$Cp) #best Cp
  Xte <- as.matrix.data.frame(dfTest)
  yte <- Xte[,ncol(Xte)]
  Xte <- Xte[,-ncol(Xte)]
  predict(ans, newx=Xte, s=iopt)$fit
}
