yhat_plus <-
function(dfTrain, dfTest, normalize=TRUE, ic=c("BIC","AIC"),
                         method=c("scad", "mc+", "lasso")){
  method <- match.arg(method)
  ic <- match.arg(ic)
  Xtr <- as.matrix.data.frame(dfTrain)
  ytr <- Xtr[,ncol(Xtr)]
  Xtr <- Xtr[,-ncol(Xtr)]
  ans <- plus(Xtr, ytr, method=method, normalize=normalize)
  n <- nrow(Xtr)
  Dev <- n*log((1-ans$r.square)/n)
  k <- ifelse(identical(ic, "AIC"), 2, log(n))
  ICs <- Dev + k*ans$dim
  lamIC <- ans$lam[which.min(ICs)]
  Xte <- as.matrix.data.frame(dfTest)
  yte <- Xte[,ncol(Xte)]
  Xte <- Xte[,-ncol(Xte)]
  capture.output(
    out<-(predict(ans, newx=Xte, lam=ans$lam)$newy)[which.min(ICs),])
  out
}
