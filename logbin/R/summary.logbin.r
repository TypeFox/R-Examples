summary.logbin <- function(object, correlation = FALSE, ...) {
  df.r <- object$df.residual
  p <- object$rank
  coef.p <- object$coefficients
  
  dispersion <- 1
  
  if(!object$boundary) {
    x <- object$x
    s <- object$prior.weights
    y <- s * object$y
    mu <- object$fitted.values
	info <- t(x) %*% apply(x,2,"*",s*mu/(1-mu))
    covmat.unscaled <- try(solve(info), silent = TRUE)
    if(!inherits(covmat.unscaled,"try-error") | all(is.nan(covmat.unscaled))) {
      covmat.scaled <- dispersion * covmat.unscaled
      var.cf <- diag(covmat.scaled)
      s.err <- sqrt(var.cf)
      tvalue <- coef.p/s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
      coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    } else {
      warning("summary.logbin: observed information matrix is singular, could not calculate covariance matrix", call. = FALSE)
      covmat.unscaled <- matrix(NaN, p, p)
      covmat.scaled <- matrix(NaN, p, p)
      coef.table <- cbind(coef.p, NaN, NaN, NaN)
    }
  }
  else {
    warning("MLE on boundary of parameter space, cannot use asymptotic covariance matrix", call. = FALSE)
    covmat.unscaled <- matrix(NaN, p, p)
    covmat.scaled <- matrix(NaN, p, p)
    coef.table <- cbind(coef.p, NaN, NaN, NaN)
  }
  
  dimnames(covmat.unscaled) <- dimnames(covmat.scaled) <- list(names(coef.p), names(coef.p))
  dimnames(coef.table) <- list(names(coef.p), c("Estimate","Std. Error","z value","Pr(>|z|)"))
  
  aliased <- rep(FALSE, p)
  names(aliased) <- names(coef.p)
  
  keep <- match(c("call", "family", "deviance", "aic", "aic.c", "df.residual",
                  "null.deviance", "df.null", "iter", "na.action"), names(object), 0L)
  ans <- c(object[keep], list(deviance.resid = residuals(object,type="deviance"),
                              coefficients = coef.table, aliased = FALSE,
                              dispersion = dispersion, df = c(p, df.r, p),
                              cov.unscaled = covmat.unscaled, cov.scaled = covmat.scaled))
  if(correlation && !any(is.nan(covmat.unscaled))) {
    dd <- sqrt(diag(covmat.unscaled))
    ans$correlation <- covmat.unscaled/outer(dd, dd)
  }
  if(inherits(object,"logbin.smooth")) ans$knots <- object$knots
  class(ans) <- c("summary.logbin", "summary.glm")
  ans
}