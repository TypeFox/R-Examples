predict.genlasso <- function(object, lambda, nlam, df, Xnew, ...) {
  # Copy the coef object and wipe out beta and u 
  co = coef(object,lambda,nlam,df)
  pr = co
  pr$beta = NULL
  pr$u = NULL
  
  if (missing(Xnew)) {
    if (is.null(object$X)) pr$fit = co$beta
    else pr$fit = object$X %*% co$beta
  }
  else {
    if (!is.matrix(Xnew) || ncol(Xnew)!=nrow(co$beta)) {
      stop(paste("Xnew must be a matrix with the appropriate number of columns.",
                 "See help pages for details."))
    }
    pr$fit = Xnew %*% co$beta
  }
  
  return(pr)
}
