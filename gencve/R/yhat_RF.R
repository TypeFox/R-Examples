yhat_RF <- function(dfTrain, dfTest) {
  Xtr <- as.matrix.data.frame(dfTrain) #assumes Gaussian!
  Xte <- as.matrix.data.frame(dfTest)
  p <- ncol(Xtr)-1
  ytr <- Xtr[, p+1]
  Xtr <- Xtr[,-(p+1)]
  Xte <- Xte[,-(p+1)]
  ans <- randomForest::randomForest(x=Xtr, y=ytr)
  predict(ans, newdata=Xte)
}

yhat_CART <- function(dfTrain, dfTest) {
  Xtr <- as.matrix.data.frame(dfTrain) #assumes Gaussian!
  Xte <- as.matrix.data.frame(dfTest)
  p <- ncol(Xtr)-1
  ytr <- Xtr[, p+1]
  Xtr <-  as.data.frame.matrix(Xtr[,-(p+1)])
  Xte <- as.data.frame.matrix(Xte[,-(p+1)])
  ans<- rpart::rpart(ytr ~., data=Xtr)
  predict(ans, newdata=Xte, type="vector")
  }
