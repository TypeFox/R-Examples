glmrob <-
function (formula, family, data, weights, subset,
	  na.action, start = NULL, offset,
          method = c("Mqle", "BY", "WBY", "MT"),
	  weights.on.x = c("none", "hat", "robCov", "covMcd"), control = NULL,
	  model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, trace.lev = 0,
	  ...)
{
    call <- match.call()
    if (is.character(family))
	family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
	family <- family()
    fami <- family$family
    if(is.null(fami))
	stop(gettextf("'%s' is not a valid family (see ?family)",
		      as.character(call[["family"]])), domain=NA)

    if (!(fami %in% c("binomial", "poisson", "Gamma", "gaussian"))) {
	stop(gettextf("Robust GLM fitting not yet implemented for family %s",
			  fami), domain=NA)
    }
    if (missing(data))
	data <- environment(formula)
    ##
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
	       names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if(identical(method, "model.frame")) return(mf)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")# "numeric" or "factor"
    if (length(dim(Y)) == 1) {
	nm <- rownames(Y)
	dim(Y) <- NULL
	if (!is.null(nm))
	    names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
	model.matrix(mt, mf, contrasts) else matrix(, NROW(Y), 0)
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0))
	stop("'weights' must be non-negative")
    if (!is.null(offset) && length(offset) != NROW(Y))
	stop(gettextf("Number of offsets is %d, should rather equal %d (number of observations)",
		      length(offset), NROW(Y)), domain=NA)
    method <- match.arg(method)
    meth. <- if(method == "WBY") "BY" else method
### FIXME: the whole 'control' should be changed to "copy"  lmrob() and lmrob.control()
## ------- --> *one* exported glmrob.control() function with 'method' and switch() inside...
## see >>> ./lmrob.MM.R

    if(is.null(control)) # -> use e.g., glmrobMqle.control()
	control <- get(paste0("glmrob", meth., ".control"))(...)
    if(missing(weights.on.x) || is.character(weights.on.x))
        weights.on.x <- match.arg(weights.on.x)
    else if(!(is.function(weights.on.x) || is.list(weights.on.x) ||
              (is.numeric(weights.on.x) && length(weights.on.x) == NROW(Y))))
        stop("'weights.on.x' must be a string, function, list or numeric n-vector")
    if(!is.null(start) && !is.numeric(start)) {
	## initialization methods
	if(!is.character(start))
	    stop("'start' must be a numeric vector, NULL, or a character string")
	start <-
	    switch(start,
		   "lmrob" =, "lmrobMM" = {
		       if(!is.null(weights))
			   warnings("weights are not yet used in computing start estimate")
		       lmrob.fit(x = X, y = family$linkinv(Y),
				 control=lmrob.control())$coefficients
		   },
		   stop("invalid 'start' string"))
    }
    fit <- switch(method,
		  "cubif" = stop("For method 'cubif', use glmRob() from package 'robust'")
		  ,
		  "Mqle" = ## --> ./glmrobMqle.R
		  glmrobMqle(X = X, y = Y, weights = weights, start = start,
			     offset = offset, family = family,
			     weights.on.x = weights.on.x, control = control,
			     intercept = attr(mt, "intercept") > 0, trace=trace.lev),
                  "BY" =, "WBY" = {
                      if(fami != "binomial")
                          stop(gettextf(
			"method='%s' is only applicable for binomial family, but family=\"\"",
                              method, fami), domain=NA)
                      ### FIXME: use glmrobBY(..) with these arguments, including 'weights'
                      glmrobBY(X=X, y=Y, weights=weights, start=start,
                               method=method, ## == "BY" / "WBY"
                               weights.on.x = weights.on.x, control = control,
                               intercept = attr(mt, "intercept") > 0,
                               trace.lev=trace.lev)
                  },
                  "MT" = {
                      glmrobMT(x=X,y=Y, weights=weights, start=start, offset = offset,
			       family=family, weights.on.x=weights.on.x, control=control,
                               intercept = attr(mt, "intercept") > 0, trace.lev=trace.lev)
                  },
		  stop("invalid 'method': ", method))
    ##-	    if (any(offset) && attr(mt, "intercept") > 0) {
    ##-		fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE],
    ##-		    y = Y, weights = weights, offset = offset,
    ##-		    control = control, intercept = TRUE)$deviance
    ##-	    }
    fit$na.action <- attr(mf, "na.action")
    if (model)
	fit$model <- mf
    if (x)
	fit$x <- X
    if (!y) ## fit$y <- NULL
	warning("setting 'y = FALSE' has no longer any effect")
    fit <- c(fit,
	     list(call = call, formula = formula, terms = mt, data = data,
		  offset = offset, control = control, method = method,
		  prior.weights = if(is.null(weights)) rep.int(1, nrow(X))
		  else weights,
		  contrasts = attr(X, "contrasts"),
		  xlevels = .getXlevels(mt, mf)))
    class(fit) <- c("glmrob", "glm")
    fit
}


summary.glmrob <- function(object, correlation=FALSE, symbolic.cor=FALSE, ...)
{
    dispersion <- object$dispersion
    if(is.null(dispersion)) dispersion <- 1
    coefs <- object$coefficients
    aliased <- is.na(coefs)# needs care; also used in print method
    if(any(aliased))
	coefs <- coefs[!aliased]
    covmat <- object$cov
    s.err <- sqrt(diag(covmat))
    zvalue <- coefs/s.err
    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind("Estimate" = coefs, "Std. Error" = s.err,
			"z value" = zvalue, "Pr(>|z|)" = pvalue)

    ans <- c(object[c("call", "terms", "family", "iter", "control", "method",
		      "residuals", "fitted.values", "w.r", "w.x")],
	     ## MM: should rather keep more from 'object' ?
	     ##	    currently, cannot even print the asympt.efficiency!
	     list(deviance=NULL, df.residual=NULL, null.deviance=NULL,
		  df.null= NULL, df= NULL, ## (because of 0 weights; hmm,...)
		  aliased = aliased,
		  coefficients = coef.table, dispersion = dispersion,
		  cov.scaled = covmat))
    if (correlation) {
	ans$correlation <- cov2cor(covmat)
	ans$symbolic.cor <- symbolic.cor
    }
    structure(ans, class = "summary.glmrob")
}

## almost a copy of vcov.glm() [if that didn't have summmary.glm() explicitly]
vcov.glmrob <- function (object, ...)
{
    so <- summary(object, corr = FALSE, ...)
    ## so$dispersion * so$cov.unscaled
    ## changed from cov.unscaled to cov.scaled
    so$cov.scaled
}


print.glmrob <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (length(coef(x))) {
	cat("Coefficients")
	if (is.character(co <- x$contrasts))
	    cat("  [contrasts: ", apply(cbind(names(co), co),
					1, paste, collapse = "="), "]")
	cat(":\n")
	print.default(format(x$coefficients, digits = digits),
		      print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nNumber of observations:", length(x$residuals),
	"\nFitted by method ", sQuote(x$method), "\n")
    invisible(x)
}

print.summary.glmrob <-
    function (x, digits = max(3, getOption("digits") - 3),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (length(cf <- coef(x))) {
	if(nsingular <- sum(x$aliased)) # glm has   df[3] - df[1]
	    cat("\nCoefficients: (", nsingular,
		" not defined because of singularities)\n", sep = "")
	else cat("\nCoefficients:\n")
	printCoefmat(cf, digits = digits, signif.stars = signif.stars,
		     na.print = "NA", ...)

	summarizeRobWeights(x$w.r * x$w.x, digits = digits,
			    header = "Robustness weights w.r * w.x:", ...)
    }
    else cat("No coefficients\n\n")

    n <- length(x$residuals)
    cat("\nNumber of observations:", n,
	"\nFitted by method", sQuote(x$method)," (in", x$iter, "iterations)\n")

    cat("\n(Dispersion parameter for ", x$family$family,
	" family taken to be ", format(x$dispersion), ")\n\n",sep = "")
    if(any(!is.null(unlist(x[c("null.deviance", "deviance")]))))
	cat(apply(cbind(paste(format(c("Null", "Residual"), justify="right"),
			      "deviance:"),
			format(unlist(x[c("null.deviance", "deviance")]),
			       digits=max(5, digits + 1)), " on",
			format(unlist(x[c("df.null", "df.residual")])),
			" degrees of freedom\n"),
		  1L, paste, collapse=" "), "\n", sep = "")
    else
	cat("No deviance values available \n")
    correl <- x$correlation
    if (!is.null(correl)) {
	p <- NCOL(correl)
	if (p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if (isTRUE(symbolic.cor)) {
		print(symnum(correl, abbr.colnames=NULL))
	    }
	    else {
		correl <- format(round(correl, 2), nsmall=2, digits=digits)
		correl[!lower.tri(correl)] <- ""
		print(correl[-1, -p, drop=FALSE], quote=FALSE)
	    }
	}
    }

    printControl(x$control, digits = digits)
    cat("\n")
    invisible(x)
}

weights.glmrob <- function(object, type = c("prior", "robustness"), ...) {
    type <- match.arg(type)
    w <- if (type == "prior") {
	## Issue warning only if called from toplevel. Otherwise the warning pop
	## up at quite unexpected places, e.g., case.names().
	if (is.null(object[["weights"]]) && identical(parent.frame(), .GlobalEnv))
	    warning("No weights defined for this object. Use type=\"robustness\" argument to get robustness weights.")
	object[["weights"]]
    } else object$w.r * object$w.x ## those also used summarizeRobWeights(x$w.r * x$w.x, ..)
    if (is.null(object$na.action)) w else naresid(object$na.action, w)
}

## Stems from a copy of residuals.glm() in
## ~/R/D/r-devel/R/src/library/stats/R/glm.R
residuals.glmrob <-
    function(object,
	     type = c("deviance", "pearson", "working", "response",
             "partial"),
	     ...)
{
    type <- match.arg(type)
    y <- object$y
    r <- object$residuals
    mu	<- object$fitted.values
    wts <- object$prior.weights # ok
    p <- length(object$coefficients)
    switch(type,
           deviance=, pearson=, response=
           if(is.null(y)) {
               mu.eta <- object$family$mu.eta
               eta <- object$linear.predictors
               ## we cannot use 'r <- ...$residuals' __ FIXME __
               stop("need non-robust working residuals for this model type")
               y <-  mu + r * mu.eta(eta)
           })
    res <- switch(type,
##		  deviance = if(object$df.residual > 0) {
		  deviance = if((nobs(object) - p) > 0) {
		      d.res <- sqrt(pmax.int((object$family$dev.resids)(y, mu, wts), 0))
		      ifelse(y > mu, d.res, -d.res)
		  } else rep.int(0, length(mu)),
		  pearson = (y-mu)*sqrt(wts)/sqrt(object$family$variance(mu)),
		  working = r,
		  response = y - mu,
		  partial = r
		  )
    if(!is.null(object$na.action))
        res <- naresid(object$na.action, res)
    if (type == "partial") ## need to avoid doing naresid() twice.
        res <- res+predict(object, type="terms")
    res
}
