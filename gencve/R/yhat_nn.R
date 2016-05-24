yhat_nn <-
function(dfTrain, dfTest, normalize=TRUE){
  nte <- nrow(dfTest)
  ntr <- nrow(dfTrain)
  p <- ncol(dfTrain)-1
  Xtr <- as.matrix.data.frame(dfTrain)
  ytr <- Xtr[,p+1]
  Xtr <- Xtr[,1:p]
  Xte <- as.matrix.data.frame(dfTest)
  Xte <- Xte[,1:p]
  if (normalize) {#rescale training. Then rescale test with same!!
    Xtr <- scale(Xtr)
    a <- attr(Xtr, "scaled:center")
    b <- attr(Xtr, "scaled:scale")
    Xte <- sweep(Xte, 2, a)
    Xte <- sweep(Xte, 2, b, FUN="/")
  }
  yHat <- numeric(nte)
  for (i in 1:nte) {
    xi <- Xte[i,]
    edist <- rowSums((Xtr-matrix(xi, byrow=TRUE, ncol=p, nrow=ntr))^2)
    yHat[i] <- ytr[which.min(edist)]
  }
yHat
}
