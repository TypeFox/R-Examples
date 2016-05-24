cov.sel.high.lasso<-function (Y, X,minscreen = 2, ...){
  alpha<-1
  family<-ifelse(length(unique(Y))>2,"gaussian","binomial")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
  }
  fitCV <- cv.glmnet(x = X, y = Y, lambda = NULL, type.measure = "deviance", family = family, alpha = alpha)
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.1se))[-1] !=0)
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
    sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2,
                     function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[,
                                                       newCut] != 0)
  }
  return(whichVariable)
}
