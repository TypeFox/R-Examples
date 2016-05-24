`deviance.mark` <- function(object, ...) object$results[['deviance']]

`confint.mark` <- function (object, parm, level = 0.95, ...) {
    cf <- object$results$beta[, 1L]
	nm <- names(cf) <- rownames(object$results$beta)
	df.residual <- object$results$n - object$results$npar
	vcv <- object$results$beta.vcv
	dimnames(vcv) <- list(nm, nm)
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    fac <- qt(a, df.residual)
    pct <- format.perc(a, 3L)
	
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
        pct))
    ses <- sqrt(diag(vcv))[parm]
    ci[] <- cf[parm] + ses %o% fac
    ci
}

`formula.mark` <- function (x, expand = TRUE, ...) {
	param <- if(is.null(x$model.parameters)) x$parameters else  x$model.parameters
	f <- lapply(param, "[[", 'formula')
	f <- f[!vapply(f, is.null, FALSE)]
	
	npty <- length(f)
	z <- vector(npty, mode = "list")
	pty <- names(f)
	
	if(expand) {
		for(i in seq_len(npty)) z[[i]] <- paste0(pty[i], "(",
				getAllTerms(f[[i]], intercept = TRUE), ")")
		res <- reformulate(gsub("((Intercept))", "(1)", unlist(z), fixed = TRUE))
	} else {
		for(i in seq_len(npty)) z[[i]] <- call(pty[i], f[[i]][[2L]])
		res <- z[[1L]]
		if(npty > 1L) for(i in seq(2L, npty)) res <- call("+", res, z[[i]])
		res <- eval(call("~", res))		
	}
	environment(res) <- environment(f[[1L]])
	res
}

