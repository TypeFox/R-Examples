yh_CART <- function(dfTr, dfTe){
  p <- ncol(dfTr)-1
  X <- dfTr[,1:p]
  y <- factor(dfTr[,p+1])
  Xte <- dfTe[,1:p]
  yte <- factor(dfTe[,p+1])
  ans<- rpart::rpart(y ~., data=X)
  yh <- predict(ans, Xte, type="class")
  YTE <- as.numeric(yte)
  YH <-  as.numeric(yh)
  if (var(YTE)*var(YH) > 0) {
    pcorr <- cor(YTE, YH)
  } else {
    pcorr <- NA
  }
  unlist(list(cost=misclassificationrate(yte, yh),
              pcorr=pcorr))
}

yh_C50 <- function(dfTr, dfTe){
  p <- ncol(dfTr)-1
  X <- dfTr[,1:p]
  y <- factor(dfTr[,p+1])
  Xte <- dfTe[,1:p]
  yte <- factor(dfTe[,p+1])
  ans<- C50::C5.0(y ~., data=X)
  yh <- C50::predict.C5.0(ans, Xte, type="class")
  YTE <- as.numeric(yte)
  YH <-  as.numeric(yh)
  if (var(YTE)*var(YH) > 0) {
    pcorr <- cor(YTE, YH)
  } else {
    pcorr <- NA
  }
  unlist(list(cost=misclassificationrate(yte, yh),
              pcorr=pcorr))
}

yh_RF <- function(dfTr, dfTe){
  p <- ncol(dfTr)-1
  X <- dfTr[,1:p]
  y <- factor(dfTr[,p+1])
  Xte <- dfTe[,1:p]
  yte <- factor(dfTe[,p+1])
  ans<- randomForest::randomForest(y ~., data=X)
  yh <- predict(ans, Xte, type="response")
  YTE <- as.numeric(yte)
  YH <-  as.numeric(yh)
  if (var(YTE)*var(YH) > 0) {
    pcorr <- cor(YTE, YH)
  } else {
    pcorr <- NA
  }
  unlist(list(cost=misclassificationrate(yte, yh),
              pcorr=pcorr))
}

yh_NB <- function(dfTr, dfTe){
  p <- ncol(dfTr)-1
  X <- dfTr[,1:p]
  y <- factor(dfTr[,p+1])
  Xte <- dfTe[,1:p]
  yte <- factor(dfTe[,p+1])
  ans<- e1071::naiveBayes(x=X, y=y)
  yh <- predict(ans, Xte, type="class")
  YTE <- as.numeric(yte)
  YH <-  as.numeric(yh)
  if (var(YTE)*var(YH) > 0) {
    pcorr <- cor(YTE, YH)
  } else {
    pcorr <- NA
  }
  unlist(list(cost=misclassificationrate(yte, yh),
              pcorr=pcorr))
}


