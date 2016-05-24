## Measures of precision and bias for ridge regression
## 

precision <- function(object, ...) {
	UseMethod("precision")
}

# DONE:  allow choice of log.det or det()^{1/p}
precision.ridge <- function(object, det.fun=c("log","root"), normalize=TRUE, ...) {
	tr <- function(x) sum(diag(x))
	maxeig <- function(x) max(eigen(x)$values)
	
	V <- object$cov
	p <- ncol(coef(object))
	det.fun <- match.arg(det.fun)
	ldet <- unlist(lapply(V, det))
	ldet <- if(det.fun == "log") log(ldet) else ldet^(1/p)
	trace <- unlist(lapply(V, tr))
	meig <- unlist(lapply(V, maxeig))	
	norm <- sqrt(rowMeans(coef(object)^2))
	if (normalize) norm <- norm / max(norm)
	data.frame(lambda=object$lambda, df=object$df, det=ldet, trace=trace, max.eig=meig, norm.beta=norm)
}

precision.lm <- function(object, det.fun=c("log","root"), normalize=TRUE, ...) {
	V <- vcov(object)
	beta <- coefficients(object)
	if (names(beta[1]) == "(Intercept)") {
		V <- V[-1, -1]
		beta <- beta[-1]
	}
#	else warning("No intercept: precision may not be sensible.")
	p <- length(beta)
	det.fun <- match.arg(det.fun)
	ldet <- det(V)
	ldet <- if(det.fun == "log") log(ldet) else ldet^(1/p)
	trace <- sum(diag(V))
	meig <- max(eigen(V)$values)
	norm <- sqrt(mean(beta^2))
	if (normalize) norm <- norm / max(norm)
	res <- list(df=length(beta), det=ldet, trace=trace, max.eig=meig, norm.beta=norm)
	unlist(res)
}
