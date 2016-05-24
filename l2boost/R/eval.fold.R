# This is a hidden function of the l2boost package.
# core cv function for generating the K folds

# @param k Index of this fold
# @param K Total number of folds to perform
# @param all.folds List of length K of observation indexes sorted into the K-fold data partitions
# @param x design matrix
# @param y response vector
# @param M Total number of l2boost iterations to perform
# @param nu l1 shrinkage parameter (0< nu <= 1)
# @param lambda l2 shrinkage parameter for elasticBoosting (lambda > 0 || lambda = NULL) 
# @param trace Show fold progress information? (T||F)
# @param type Which l2boost algorithm? 
# @param ... extra arguments passed into the \code{\link{l2boost}} method.
# 
# @seealso \code{\link{l2boost}}
eval.fold <- function(k, K, all.folds, x, y, M, nu, lambda, trace, type, ...) { 
  if (trace) {
    if (k <= K) {
      cat("\t K-fold:", k, "\n")
    }
    else {
      cat("\t final analysis (full-data)\n")
    }
  }
  omit <- all.folds[[k]]
  
  fit <- l2boost(x = as.matrix(x[-omit,, drop = FALSE]), y = y[-omit],
                 M = M, nu = nu, type = type, lambda = lambda, ...=...)
  
  #print(fit)
  if (k <= K) {
    yhat.path <- predict.l2boost(fit, xnew = x[omit, , drop=FALSE])$yhat.path
    mse <- sapply(1:length(yhat.path), function(m) {
      mean((yhat.path[[m]] - y[omit])^2, na.rm = TRUE)})
    return(list(obj = fit, mse = mse)) 
  }
  else {
    return(list(obj = fit))    
  }
}
