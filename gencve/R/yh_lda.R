yh_lda <- function(dfTr, dfTe){
  p <- ncol(dfTr)-1
  X <- dfTr[,1:p]
  y <- dfTr[,p+1]
  Xte <- dfTe[,1:p]
  yte <- dfTe[,p+1]
  ans <- MASS::lda(x=X, grouping=y)
  yh <- predict(ans, newdata=Xte)$class
  unlist(list(cost=misclassificationrate(yte, yh),
              pcorr=cor(as.numeric(yte), as.numeric(yh))))
}
