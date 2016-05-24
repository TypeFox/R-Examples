predict.hqreg <- function(object, X, lambda, type=c("response","coefficients","nvars"), exact = FALSE, ...) {
  type=match.arg(type)
  if (missing(X) && type == "response") stop("Need to supply 'X'")
  beta <- coef.hqreg(object, lambda = lambda, exact = exact)
  if (type == "coefficients") return(beta)
  if (is.matrix(beta)) {
    b0 <- beta[1,]
    b <- beta[-1,]
  } else {
    b0 <- beta[1]
    b <- beta[-1]
  }
  if (type == "nvars") {
    if (is.matrix(beta)) return(apply(b!=0, 2, sum))
    else return(sum(b!=0))
  }
  if (type == "response") return(sweep(X %*% b, 2, b0, "+"))
}

coef.hqreg <- function(object, lambda, exact = FALSE, ...) {
  if (missing(lambda)) {
    beta <- object$beta
  } else if (exact) {
    # augment the lambda sequence with the new values, and refit the model
    ls <- object$lambda
    ind <- match(lambda,ls,0)
    if (any(ind == 0)) {
      ls <- unique(rev(sort(c(lambda,ls))))
      object <- update(object,lambda=ls)
    }
    beta <- object$beta[, ls == lambda]
  } else {
    # use linear interpolation to estimate coefficients for supplied lambda
    ls <- object$lambda
    lambda[lambda>max(ls)] <- max(ls)
    lambda[lambda<min(ls)] <- min(ls)
    ind <- approx(ls,seq(ls),lambda)$y
    left <- floor(ind)
    right <- ceiling(ind)
    weight <- ind %% 1
    beta <- (1-weight)*object$beta[,left] + weight*object$beta[,right]
    if (length(lambda) > 1) colnames(beta) <- round(lambda,4)
  }
  return(beta)
}
