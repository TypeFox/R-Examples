# last modified 2012-07-29 by J. Fox

raw.moments <- function(...){
	.Deprecated("rawMoments", package="sem")
	rawMoments(...)
}

rawMoments <- function(object, ...) UseMethod("rawMoments")

rawMoments.formula <- function(formula, data, subset, na.action, 
		contrasts = NULL, ...) {
	if (missing(na.action))
		na.action <- options()$na.action
	m <- match.call(expand.dots = FALSE)
	if (is.matrix(eval(m$data, sys.frame(sys.parent()))))
		m$data <- as.data.frame(data)
	m$instruments <- m$contrasts <- NULL
	m[[1]] <- as.name("model.frame")
	mf <- eval(m, sys.frame(sys.parent()))
	response <- attr(attr(mf, "terms"), "response")
	if (response) stop("formula cannot have a response")
	na.act <- attr(mf, "na.action")
	X <- model.matrix(formula, data = mf, contrasts)
	nms <- colnames(X)
	if ("(Intercept)" %in% nms) colnames(X)[nms == "(Intercept)"] <- "Intercept"
	rawMoments(X)
}

rawMoments.default <- function(object, na.rm=FALSE, ...){
	object <- as.matrix(object)
	if (na.rm) object <- na.omit(object)
	N <- nrow(object)
	result <- crossprod(object, object)/N
	attr(result, "N") <- N
	class(result) <- "rawmoments"
	result
}

print.rawmoments <- function(x, ...){
	xx <- unclass(x)
	attr(xx, "N") <- NULL
	cat("\nRaw Moments\n")
	print(xx, ...)
	cat("\nN = ", attr(x, "N"), "\n")
	invisible(x)
}

cov2raw <- function(cov, mean, N, sd){
	if (all(1 == diag(cov)) && !missing(sd))
		cov <- cov * outer(sd, sd)
	raw <- ((N - 1)*cov + N*outer(mean, mean))/N
	colnames(raw) <- rownames(raw) <- rownames(cov)
	raw <- rbind(c(1, mean), cbind(mean, raw))
	if (!("Intercept" %in% rownames(raw)))
		colnames(raw)[1] <- rownames(raw)[1] <- "Intercept"
	attr(raw, "N") <- N
	class(raw) <- "rawmoments"
	raw
} 
