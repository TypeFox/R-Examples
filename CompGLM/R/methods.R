#' @export
print.Comp <- function(x, ...) {	
	cat("\nCall:\n")
	print(x$call)
	cat("\nBeta:\n")
	print(x$beta)
	cat("\nZeta:\n")
	print(x$zeta)
	cat("\nDegrees of Freedom:", x$df.null,"Total (i.e. Null);", x$df.residual , "Residual\n")
	cat("AIC:", AIC(x),"\n")
	cat("Log-Likelihood:", x$logLik, "\n")
	invisible(x)
}

#' @export
summary.Comp <- function(object, ...) {	
	vcov <- vcov(object)
	se <- sqrt(diag(vcov))
	
	betaCoefs <- object$beta
	zetaCoefs <- object$zeta
	
	seBeta <- se[seq_along(betaCoefs)]
	seZeta <- se[-seq_along(betaCoefs)]
	
	tValuesBeta <- betaCoefs / seBeta
	tValuesZeta <- zetaCoefs / seZeta
	
	coefTableBeta <- cbind(Estimate = betaCoefs, Std.Error = seBeta, t.value = tValuesBeta, 
			p.value = 2.0 * pt(-abs(tValuesBeta), df = object$df.residual))
	coefTableZeta <- cbind(Estimate = zetaCoefs, Std.Error = seZeta, t.value = tValuesZeta, 
			p.value = 2.0 * pt(-abs(tValuesZeta), df = object$df.residual))
	
	res <- list(call = object$call, beta = coefTableBeta, zeta = coefTableZeta, 
			AIC = AIC(object), logLik = object$logLik)
	
	class(res) <- "summary.Comp"
	return(res)
}

#' @export
print.summary.Comp <- function(x, ...) {	
	cat("\nCall:\n")
	print(x$call)
	cat("\nBeta:\n")
	printCoefmat(x$beta, P.values = TRUE, has.Pvalue = TRUE, 
			signif.stars = TRUE, signif.legend = FALSE)
	cat("\nZeta:\n")
	printCoefmat(x$zeta, P.values = TRUE, has.Pvalue = TRUE, 
			signif.stars = TRUE, signif.legend = TRUE)
	cat("\nAIC:", x$AIC, "\n")
	cat("Log-Likelihood:", x$logLik, "\n")
	invisible(x)
}

#' @export
vcov.Comp <- function(object, ...) {	
	vcov <- solve(object$hessian)
	return(vcov)
}

#' @export
logLik.Comp <- function(object, ...) {	
	ans <- object$logLik
	df <- length(c(object$beta, object$zeta))
	attr(ans, "df") <- df
	attr(ans, "nobs") <- nobs(object)
	class(ans) <- "logLik"
	return(ans)
}

#' @export
coef.Comp <- function(object, ...) {
	return(list(beta = object$beta, zeta = object$zeta))
}
 
#' @export
extractAIC.Comp <- function(fit, scale, k, ...) {
	edf <- length(fit$coefficients)
	AIC <- -2.0 * fit$loglikelihood + 2.0 * edf
	ans <- c(edf, AIC)
	return(ans)
}

#' @export
nobs.Comp <- function(object, ...) {
	return(object$nobs)
}

#' @export
predict.Comp <- function(object, newdata, ...) {
	
	if (missing(newdata)) {
		newdata <- object$data
	}
	
	modelTerms <- delete.response(object$terms)
	modelFrame <- model.frame(modelTerms, newdata)
	
	beta <- object$beta
	lamOffset <- object$lamOffset
	xLam <- model.matrix(modelTerms, modelFrame)
	lam <- exp(xLam %*% beta + lamOffset)
	
	if (is.null(object$nuModelTerms)) {
		xNu <- matrix(1.0, nobs(object), 1L)
	}	else {
		nuModelTerms <- object$nuModelTerms
		nuModelFrame <- model.frame(nuModelTerms, newdata)
		xNu <- model.matrix(nuModelTerms, nuModelFrame)
	}
	zeta <- object$zeta
	nuOffset <- object$nuOffset
	nu <- exp(xNu %*% zeta + nuOffset)
	
	out <- matrix(NA_real_, nobs(object), 2L)
	out[ , 1L] <- lam
	out[ , 2L] <- nu
	colnames(out) <- c("lambda", "nu")
	
	return(out)
}
