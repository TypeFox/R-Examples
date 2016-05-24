softthresh <- function(object, lambda, gamma) {
  if (!is.numeric(gamma) || gamma<0) {
    stop("gamma must be nonnegative.")
  }
  if (!any(class(object)=="fusedlasso") || !is.null(object$X) ||
      (!is.null(object$gamma) && object$gamma!=0)) {
    warning(paste("Soft-thresholding only gives a valid primal solution when applied",
            "to a fused lasso problem with pure fusion (gamma=0) and identity predictor matrix X."))
  }

  # If no lambdas were passed, take the knots in the path
  if (missing(lambda)) lambda = object$lambda

  beta = coef(object,lambda=lambda)$beta
  lams = matrix(gamma*lambda,nrow(beta),ncol(beta),byrow=TRUE)
  beta = sign(beta)*pmax(abs(beta)-lams,0)

  return(beta)
}
