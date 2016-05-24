yhat_step <-
function(dfTrain, dfTest, ic=c("BIC","AIC")) {
    ansFull <- lm(y~., data=dfTrain)
    ic <- match.arg(ic)
    k <- ifelse(identical(ic, "AIC"), 2, log(nrow(dfTrain)))
    junk <- capture.output(ansStep <- step(ansFull, k=k))
    predict(ansStep, newdata=dfTest)
  }
