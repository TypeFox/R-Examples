yh_kNN <- function(dfTr, dfTe,
                   method=c("LOOCV", "MLE", "NN"), k=1) {
  method<- match.arg(method)
  p <- ncol(dfTr)-1
  X <- dfTr[,1:p]
  y <- dfTr[,p+1]
  Xte <- dfTe[,1:p]
  yte <- dfTe[,p+1]
  kOpt <- switch(method,
            LOOCV = kNN_LOOCV(X, y),
            MLE = kNN_MLE(X, y),
            NN = k)
  yh <- class::knn(X, Xte, y, k=kOpt)
  unlist(list(cost=misclassificationrate(yte, yh),
              pcorr=cor(as.numeric(yte), as.numeric(yh))))
}

yh_NN <- function(dfTr, dfTe) {
  p <- ncol(dfTr)-1
  X <- dfTr[,1:p]
  y <- dfTr[,p+1]
  Xte <- dfTe[,1:p]
  yte <- dfTe[,p+1]
  yh <- class::knn1(train=X, test=Xte, cl=y)
  unlist(list(cost=misclassificationrate(yte, yh),
              pcorr=cor(as.numeric(yte), as.numeric(yh))))
}


