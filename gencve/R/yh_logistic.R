yh_logistic <- function(dfTr, dfTe, alpha=NULL){
  fixupFactor <- function(u) factor(as.character(u)) #needed, eg: churn etc
  p <- ncol(dfTr)-1
  X <- dfTr[,1:p]
  y <- fixupFactor(dfTr[,p+1])
  stopifnot(length(levels(y))==2)
  Xte <- dfTe[,1:p]
  yte <-  fixupFactor(dfTe[,p+1])
  if (is.null(alpha)) {#glm
    ans<- glm(y ~., data=X, family="binomial")
    prob <- ifelse(predict.glm(ans, newdata=Xte, type="response")<0.5, 1, 2)
    yh <- levels(y)[prob]
  } else {#elastic net
    stopifnot(alpha<=1 || alpha>=0)
    X <- as.matrix.data.frame(X)
    Xte <- as.matrix.data.frame(Xte)
    ans<- glmnet(x=X, y=y, family="binomial", alpha=alpha)
    ans_cv <- cv.glmnet(X, y=y, family="binomial", alpha=alpha)
    lambdaHat <- ans_cv$lambda.1se
    yh <- predict(ans, newx=Xte, type="class", s=lambdaHat)[,1]
  }
  yh <- factor(yh)
  unlist(list(cost=misclassificationrate(yte, yh),
                     pcor=cor(as.numeric(yte), as.numeric(yh))))
}


