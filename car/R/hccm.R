#-------------------------------------------------------------------------------
# Revision history:
# 2009-09-16: optionally allow models with aliased coefficients J. Fox
# 2012-04-04: modified to allow weighted linear models. J. Fox
#-------------------------------------------------------------------------------

# Heteroscedasticity-corrected standard errors (White adjustment) (J. Fox)

hccm <- function(model, ...){
	#last modified 12 Dec 2000 by J. Fox
	UseMethod("hccm")
}

#hccm.lm <- function (model, type = c("hc3", "hc0", "hc1", "hc2", "hc4"), 
#		singular.ok=TRUE, ...) {
#	if (!is.null(weights(model))) 
#		stop("requires unweighted lm")
#	type <- match.arg(type)
#	if (any(aliased <- is.na(coef(model))) && !singular.ok)
#		stop("there are aliased coefficients in the model")
#	sumry <- summary(model, corr = FALSE)
#	s2 <- sumry$sigma^2
#	V <- sumry$cov.unscaled
#	if (type == FALSE) 
#		return(s2 * V)
#	e <- na.omit(residuals(model))
#	X <- model.matrix(model)[, !aliased]
#	df.res <- df.residual(model)
#	n <- length(e)
#	h <- hat(X)
#	p <- ncol(X)
#	factor <- switch(type, hc0 = 1, hc1 = df.res/n, hc2 = 1 - 
#					h, hc3 = (1 - h)^2, hc4 = (1 - h)^pmin(4, n * h/p))
#	V %*% t(X) %*% apply(X, 2, "*", (e^2)/factor) %*% V
#}

hccm.lm <-function (model, type = c("hc3", "hc0", "hc1", "hc2", "hc4"), 
		singular.ok = TRUE, ...) {
	e <- na.omit(residuals(model))
	removed <- attr(e, "na.action")
	wts <- if (is.null(weights(model))) 1 
			else weights(model)
	type <- match.arg(type)
	if (any(aliased <- is.na(coef(model))) && !singular.ok) 
		stop("there are aliased coefficients in the model")
	sumry <- summary(model, corr = FALSE)
	s2 <- sumry$sigma^2
	V <- sumry$cov.unscaled
	if (type == FALSE) 
		return(s2 * V)
	h <- hatvalues(model)
	if (!is.null(removed)){
		wts <- wts[-removed]
		h <- h[-removed]
	}
	X <- model.matrix(model)[, !aliased]
	df.res <- df.residual(model)
	n <- length(e)
	e <- wts*e
	p <- ncol(X)
	factor <- switch(type, hc0 = 1, hc1 = df.res/n, 
			hc2 = 1 - h, hc3 = (1 - h)^2, hc4 = (1 - h)^pmin(4, n * h/p))
	V %*% t(X) %*% apply(X, 2, "*", (e^2)/factor) %*% V
}

hccm.default<-function(model, ...){
	#last modified 12 Dec 2000 by J. Fox
	stop("requires an lm object")
}
