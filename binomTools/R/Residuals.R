exact.deletion <- function(object) {

  if(is.null(object$model)) 
    stop("glm object has to be fitted with \'model=TRUE\'")

  start <- coef(object) # start values to be used in model fitting
  fit <- object$fitted
  devResidAbs <- abs(deviance(object))
  mf <- model.frame(object)
  N <- nrow(mf)

  if(!is.null(model.weights(mf))) mw <- model.weights(mf)
  else mw <- rep(1, nrow(object$model))

  resp <- model.response(mf)
  if(NCOL(resp) == 1L) {
    n <- rep(1, nrow(mf))
    resp <- cbind(resp, n-resp)
    y <- resp
  }
  else y <- resp[,1] / rowSums(resp)

  object$model[,1] <- resp*mw		### new response including the weights
  newCall <- update(object, subset=-i, start=start, evaluate=F)

  residLik <- rep(0, N)
  for(i in 1:N) {
    objLeave1out <- eval(newCall)
    devResidLeave1out <- deviance(objLeave1out)
    residLik[i] <- sign(y[i]-fit[i]) * sqrt(devResidAbs-devResidLeave1out) 
  }

  names(residLik) <- 1:N
  residLik
}


Residuals <- function(object, type=c("approx.deletion",
                      "exact.deletion", "standard.deviance", "standard.pearson", 
                      "deviance", "pearson", "working", "response", "partial"))
{  
  type <- match.arg(type)

  res <- switch(type, 
    approx.deletion = rstudent(object), exact.deletion =
    exact.deletion(object), standard.deviance = rstandard(object),
    standard.pearson = rstandard(object, type="pearson"), deviance =
    residuals(object), pearson = residuals(object, type="pearson"),
    working = residuals(object, type="working"), response =
    residuals(object, type="response"), partial = residuals(object, type="partial"))

  return(res)
}

