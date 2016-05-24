# Code based on 'compar.gee' from package 'ape'
## Comparative Analysis with GEEs
## compar.gee.R (2011-06-14)
## Copyright 2002-2010 Emmanuel Paradis
## https://svn.mpl.ird.fr/ape/dev/ape/R/compar.gee.R

##=============================================================================
## quasiLik 
##=============================================================================

`quasiLik` <- function (object, ...) UseMethod("quasiLik")

.qlik <- function(y, mu, fam) {
	ret <- switch(fam,
		   gaussian = -sum((y - mu)^2)/2,
		   binomial = sum(y * log(mu/(1 - mu)) + log(1 - mu)),
		   #binomial.sqvar = sum(((2 * y - 1) * log(mu /(1 - mu))) - (y / mu) - ((1 - y)/(1 - mu))),
		   poisson = sum(y * log(mu) - mu),
		   Gamma = -sum(y/mu + log(mu)),
		   inverse.gaussian = sum(-y/(2 * mu^2) + 1/mu),
		   cry(, "do not know how to calculate quasi-likelihood for family '%s'",
				fam))
	ret
}

`print.quasiLik` <- function (x, digits = getOption("digits"), ...) {
    cat("'quasi Lik.' ", paste(format(c(x), digits = digits), collapse = ", "), 
        "\n", sep = "")
    invisible(x)
}

`quasiLik.geeglm` <-
`quasiLik.gee` <-
function(object, ...) {
	ret <- .qlik(object$y, object$fitted.values, family(object)$family)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(object$y)
	class(ret) <- "quasiLik"
	ret
}

`quasiLik.yagsResult` <- function(object, ...) {
	mu <- object@fitted.values
	ret <- .qlik(mu + object@residuals, mu, family(object)$family)
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(mu)
	class(ret) <- "quasiLik"
	ret
}

`quasiLik.geem` <-
function(object, ...) {
	fam <- family(object)
	ret <- .qlik(object$y, fitted(object), if(inherits(fam, "family"))
				 fam$family else "custom")
	attr(ret, "df") <- NA
	attr(ret, "nobs") <- length(object$y)
	class(ret) <- "quasiLik"
	ret
}

##=============================================================================
## QIC 
##=============================================================================

.qic2 <- function(y, mu, vbeta, mui, vbeta.naiv.i, fam, typeR = FALSE) {
	ql <- if(typeR) .qlik(y, mu, fam) else .qlik(y, mui, fam)
	# XXX: should be typeR = TRUE for QICu???
	n <- length(y)
	# yags/yags.cc: p140 of Hardin and Hilbe
	if(fam == "gaussian") ql <- (n * log(-2 * ql / n)) / -2
	AIinv <- solve(vbeta.naiv.i)
	tr <- sum(matmult(AIinv, vbeta, diag.only = TRUE)) ## TODO: use matmult
	## tr <- sum(diag(AIinv %*% vbeta)) ## TODO: use matmult
	px <- length(mu)
	## When all modelling specifications in GEE are correct tr = px.
	c(2 * (c(QIC = tr, QICu = px) - ql), n = n)
}

`getQIC` <- 
function(x, typeR = FALSE) UseMethod("getQIC")
	
`getQIC.default` <-
function(x, typeR = FALSE) .NotYetImplemented()

`getQIC.coxph` <- function(x, ...) {
	warning("QIC for 'coxph' is experimental")
	naive.var <- x[[ if (is.null(x$naive.var)) "var" else "naive.var" ]]
	# tr <- sum(diag(solve(naive.var) %*% x$var))
	tr <- sum(matmultdiag(solve(naive.var), x$var))
	ll <- x$loglik[2L]
	px <- x$n
	c(2 * (c(QIC = tr, QICu = px) - ll), n = px)
}

`getQIC.gee` <- 
function(x, typeR = FALSE) {
	if(x$model$corstr != "Independent")
		utils::capture.output(suppressMessages(xi <- update(x, corstr = "independence",
		silent = TRUE))) else
		xi <- x
	
	.qic2(x$y, x$fitted.values, x$robust.variance, 
		  xi$fitted.values, xi$naive.variance, family(x)$family,
		  typeR = typeR)
}

`getQIC.geeglm` <- 
function(x, typeR = FALSE) {
	xi <- if(x$corstr != "independence")
		update(x, corstr = "independence") else x
	.qic2(x$y, x$fitted.values, x$geese$vbeta, 
		  xi$fitted.values, xi$geese$vbeta.naiv, family(x)$family,
		  typeR = typeR)
}

`getQIC.yagsResult` <- 
function(x, typeR = FALSE) {
	xi <- if(x@corstruct.tag != "independence")
		update(x, corstruct = "independence") else x
	.qic2(x@fitted.values + x@residuals, x@fitted.values, x@robust.parmvar, 
		  xi@fitted.values, xi@naive.parmvar, family(x)$family,
		  typeR = typeR)
}

`getQIC.geem` <-
function(x, typeR = FALSE) {
	fam <- family(x)
	xi <- if(x$corr != "independence")
		update(x, corstr = "independence") else x
    .qic2(x$y, fitted(x), x$var, fitted(xi), xi$naiv.var,
		if(inherits(fam, "family")) fam$family else "custom",
        typeR = typeR)
}

`QIC` <- function (object, ..., typeR = FALSE) {
	if (!missing(...)) {
		res <- sapply(list(object, ...), getQIC, typeR = typeR)
		val <- as.data.frame(t(res[1L,, drop = FALSE]))
		colnames(val) <- c("QIC")
		Call <- match.call()
		Call$typeR <- NULL
		row.names(val) <- as.character(Call[-1L])
		val
	} else getQIC(object, typeR = typeR)[1L]
}

`QICu` <- function (object, ..., typeR = FALSE) {
	if (!missing(...)) {
		res <- sapply(list(object, ...), getQIC, typeR = typeR)
		val <- as.data.frame(t(res[2L,, drop = FALSE]))
		colnames(val) <- "QICu"
		Call <- match.call()
		Call$typeR <- NULL
		row.names(val) <- as.character(Call[-1L])
		val
	} else getQIC(object, typeR = typeR)[2L]
}
