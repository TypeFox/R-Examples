summary.addreg <- function(object, correlation = FALSE, ...) {
  df.r <- object$df.residual
  p <- object$rank - as.numeric(substr(object$family$family,1,7) == "negbin1")
  coef.p <- object$coefficients
  dispersion <- 1
  
  if(!object$boundary) {
    if(object$family$family == "poisson") {
      x <- object$x
      s <- object$standard
      y <- object$y
      eta <- object$linear.predictors
      info <- t(x) %*% apply(x,2,"*",s/eta)
      covmat.unscaled <- try(solve(info), silent = TRUE)
    } else if(object$family$family == "binomial") {
	  x <- object$x
	  s <- object$prior.weights
	  y <- s * object$y
	  eta <- object$linear.predictors
	  info1 <- t(x) %*% apply(x,2,"*",s/eta)
	  info2 <- t(x) %*% apply(x,2,"*",s/(1-eta))
	  info <- info1 + info2
	  covmat.unscaled <- try(solve(info), silent = TRUE)
	} else if(substr(object$family$family,1,7) == "negbin1") {
	  x <- object$x
	  s <- object$standard
	  y <- object$y
	  mu <- object$fitted.values
	  phi <- object$scale - 1
	  r <- mu / phi
	  dgd <- tgd <- rep(0, length(y))
	  dgd[y > 0] <- mapply(function(r, y) sum(1/(r + seq_len(y) - 1)), r = r[y > 0], y = y[y > 0])
	  tgd[y > 0] <- mapply(function(r, y) -sum(1/(r + seq_len(y) - 1)^2), r = r[y > 0], y = y[y > 0])
	  info1 <- -t(apply(x, 2, "*", s^2 * tgd)) %*% x / phi^2
	  info2 <- t(x) %*% (s/phi^2 * (r * tgd + dgd + phi/(phi+1) - log(phi+1)))
	  info3 <- -sum(r/phi * (2/phi * (dgd + phi/(phi+1) - log(phi+1)) + r/phi * tgd + phi/(phi+1)^2) - y*(2*phi+1)/(phi*(phi+1))^2)
	  info <- rbind(cbind(info1, info2), c(info2, info3))
      covmat.full <- try(solve(info), silent = TRUE)
	  if(!inherits(covmat.full,"try-error") | all(is.nan(covmat.full))) {
		covmat.unscaled <- covmat.full[(1:p),(1:p), drop = FALSE]
		var.phi <- covmat.full[(p+1),(p+1)]
	  } else {
		covmat.unscaled <- matrix(NaN, p, p)
		var.phi <- NaN
	  }
    }
    if(!inherits(covmat.unscaled,"try-error") | all(is.nan(covmat.unscaled))) {
      covmat.scaled <- dispersion * covmat.unscaled
      var.cf <- diag(covmat.scaled)
      s.err <- sqrt(var.cf)
      tvalue <- coef.p/s.err
      pvalue <- 2 * pnorm(-abs(tvalue))
      coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    } else {
      warning("summary.addreg: information matrix is singular, could not calculate covariance matrix", call. = FALSE)
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
	if(substr(object$family$family,1,7) == "negbin1") {
        phi <- object$scale - 1
        var.phi <- NaN
    }
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
  if(substr(object$family$family,1,7) == "negbin1") {
	ans$phi <- phi
	ans$var.phi <- var.phi
  }
  if(inherits(object,"addreg.smooth")) ans$knots <- object$knots
  class(ans) <- c("summary.addreg", "summary.glm")
  ans
}