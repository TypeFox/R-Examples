cv.hqreg <- function(X, y, ..., nfolds=10, fold.id, type.measure = c("deviance", "mse", "mae"), 
                     seed, trace=FALSE) {
  type.measure = match.arg(type.measure)
  if (!missing(seed)) set.seed(seed)
  fit <- hqreg(X, y, ...)
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$gamma <- fit$gamma
  loss.args <- list(method=fit$method, gamma=fit$gamma, tau=fit$tau, type.measure = type.measure)
  E <- matrix(NA, nrow=length(y), ncol=length(cv.args$lambda))
  n <- length(y)
  if(missing(fold.id)) fold.id <- ceiling(sample(1:n)/n*nfolds)
  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #",i,sep="","\n")
    cv.args$X <- X[fold.id!=i,]
    cv.args$y <- y[fold.id!=i]
    X2 <- X[fold.id==i,]
    y2 <- y[fold.id==i]
    fit.i <- do.call("hqreg", cv.args) # ensure the cross validation uses the same gamma for huber loss
    yhat <- predict(fit.i, X2)
    E[fold.id==i, 1:ncol(yhat)] <- loss.hqreg(y2, yhat, loss.args)
  }

  ## Eliminate saturated lambda values
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[,ind]
  lambda <- cv.args$lambda[ind]

  ## Results
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, sd) / sqrt(n)
  index.min <- which.min(cve)
  # adjust the selection using 1-SD method
  index.1se <- min(which(cve < cve[index.min]+cvse[index.min]))
  val <- list(cve=cve, cvse=cvse, lambda=lambda, fit=fit, lambda.1se = lambda[index.1se], lambda.min=lambda[index.min])
  structure(val, class="cv.hqreg")
}
