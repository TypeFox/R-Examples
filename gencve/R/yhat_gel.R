yhat_gel <- function(dfTrain, dfTest, alpha=1) {
#Gaussian, elastic net
  n <- nrow(dfTrain)
  n <- nrow(dfTrain)
  Xtr <- as.matrix.data.frame(dfTrain) #assumes Gaussian!
  Xte <- as.matrix.data.frame(dfTest)
  p <- ncol(Xtr)-1
  ytr <- Xtr[, p+1]
  Xtr <- Xtr[,-(p+1)]
  Xte <- Xte[,-(p+1)]
  ans<- glmnet(x=Xtr, y=ytr, alpha=alpha)
  ans_cv <- cv.glmnet(Xtr, y=ytr, alpha=alpha)
  lambdaHat <- ans_cv$lambda.1se
  predict(ans, newx=Xte, s=lambdaHat)[,1]
}
