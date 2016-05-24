yh_svm <- function(dfTr, dfTe){
  p <- ncol(dfTr)-1
  X <- dfTr[,1:p]
  y <- factor(dfTr[,p+1])
  Xte <- dfTe[,1:p]
  yte <- factor(dfTe[,p+1])
  ans <- svm(x=X, y=y)
  yh <- predict(ans, newdata=Xte)#predict.svm not exported!!
  unlist(list(cost=misclassificationrate(yte, yh),
              pcorr=cor(as.numeric(yte), as.numeric(yh))))
}
