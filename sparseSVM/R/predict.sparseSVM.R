predict.sparseSVM <- function(object, X, lambda, type=c("class","coefficients","nvars"), exact = FALSE, ...) {
  type=match.arg(type)
  if (missing(X) && type == "class") stop("Need to supply 'X'")
  weights <- coef.sparseSVM(object, lambda = lambda, exact = exact)
  if (type == "coefficients") return(weights)
  if (is.matrix(weights)) {
    b <- weights[1,]
    w <- weights[-1,]
  } else {
    b <- weights[1]
    w <- weights[-1]
  }
  if (type == "nvars") {
    if (is.matrix(w)) return(apply(w!=0, 2, sum))
    else return(sum(w!=0))
  }
  if (type == "class") {
    v <- sweep(X %*% w, 2, b, "+")
    return(ifelse(v > 0, object$levels[1], object$levels[2]))
  }
}

coef.sparseSVM <- function(object, lambda, exact = FALSE, ...) {
  if (missing(lambda)) {
    weights <- object$weights
  } else if (exact) {
    # augment the lambda sequence with the new values, and refit the model
    ls <- object$lambda
    ind <- match(lambda,ls,0)
    if (any(ind == 0)) {
      ls <- unique(rev(sort(c(lambda,ls))))
      object <- update(object,lambda=ls)
    }
    weights <- object$weights[, ls == lambda]
  } else {
    # use linear interpolation to estimate coefficients for supplied lambda
    ls <- object$lambda
    lambda[lambda>max(ls)] <- max(ls)
    lambda[lambda<min(ls)] <- min(ls)
    ind <- approx(ls,seq(ls),lambda)$y
    left <- floor(ind)
    right <- ceiling(ind)
    delta <- ind %% 1
    weights <- (1-delta)*object$weights[,left] + delta*object$weights[,right]
    if (length(lambda) > 1) colnames(weights) <- round(lambda,4)
  }
  return(weights)
}
