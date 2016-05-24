powell <- function(par, fn, control = powell.control(),
                   check.hessian = TRUE, ...) {
  dots <- list(...)
  expr <- Quote(fn(..par..))
  if(length(dots)) {
    expr[names(dots)] <- dots
    expr <- as.call(expr)
  }
  env <- new.env(list(fn = fn))
  control <- do.call("powell.control", control)
  any.na <- sapply(lapply(control, is.na), sum)
  if(any(is.na(par)) || any(any.na > 0))
    stop("missing values not allowed in call")
  n <- length(par)
  if(length(control$parscale) != n) {
    if(length(control$parscale) == 1)
      control$parscale <- rep(control$parscale, n)
    else
      stop("control$parscale must be determined for all values in \"par\"")
  }
  ret <- .Call("R_UObyQA", as.double(par), expr, control, env, PACKAGE = "powell")
  names(ret) <- c("par", "value", "counts", "hessian")
  npar <- length(ret$par)
  ret$hessian <- matrix(ret$hessian, npar, npar)
  if(!is.null(parnm <- names(par))) {
    names(ret$par) <- parnm
    dimnames(ret$hessian) <- list(parnm, parnm)
  }
  if(check.hessian) {
    ret$eigen.hessian <- eigen(ret$hessian, symmetric = TRUE)
    ev <- ret$eigen.hessian$values
    if(!all(ev >= -.Machine$double.eps^0.5 * abs(ev[1])))
      warning(ret$message <- "Hessian is not positive definite.")
  }
  names(ret$counts) <- "function"
  ret$convergence <- ifelse(control$maxit == ret$counts, 1, 0)
  ret$control <- control
  ret$call <- match.call()
  ret
}
