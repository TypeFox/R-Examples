
### The first part of lmrob()  much cut'n'paste from lm() - on purpose!
lmrob <-
    function(formula, data, subset, weights, na.action, method = 'MM',
	     model = TRUE, x = !control$compute.rd, y = FALSE,
	     singular.ok = TRUE, contrasts = NULL, offset = NULL,
	     control = NULL, init = NULL, ...)
{
    ## to avoid problems with setting argument
    ## call lmrob.control here either with or without method arg.
    if (miss.ctrl <- missing(control))
	control <- if (missing(method))
	    lmrob.control(...) else lmrob.control(method = method, ...)
    else ## check dots
	chk.s(...)
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
	       names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms") # allow model.frame to update it
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if(!is.null(w) && !is.numeric(w))
	stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if(!is.null(offset) && length(offset) != NROW(y))
	stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
		      length(offset), NROW(y)), domain = NA)
    if (!miss.ctrl && !missing(method) && method != control$method) {
	warning("Methods argument set by method is different from method in control\n",
		"Using the former, method = ", method)
	control$method <- method
    }

    if (is.empty.model(mt)) {
	x <- NULL
	singular.fit <- FALSE ## to avoid problems below
	z <- list(coefficients = if (is.matrix(y)) matrix(,0,3) else numeric(0),
		  residuals = y, scale = NA, fitted.values = 0 * y,
		  cov = matrix(,0,0), weights = w, rank = 0,
		  df.residual = NROW(y), converged = TRUE, iter = 0)
	if(!is.null(offset)) {
	    z$fitted.values <- offset
	    z$residuals <- y - offset
	    z$offset <- offset
	}
    }
    else {
	x <- model.matrix(mt, mf, contrasts)
	contrasts <- attr(x, "contrasts")
	assign <- attr(x, "assign")
	p <- ncol(x)
	if(!is.null(offset))
	    y <- y - offset
	if (!is.null(w)) {
	    ## checks and code copied/modified from lm.wfit
	    ny <- NCOL(y)
	    n <- nrow(x)
	    if (NROW(y) != n | length(w) != n)
		stop("incompatible dimensions")
	    if (any(w < 0 | is.na(w)))
		stop("missing or negative weights not allowed")
	    zero.weights <- any(w == 0)
	    if (zero.weights) {
		save.r <- y
		save.w <- w
		save.f <- y
		ok <- w != 0
		nok <- !ok
		w <- w[ok]
		x0 <- x[nok, , drop = FALSE]
		x  <- x[ ok, , drop = FALSE]
		n <- nrow(x)
		y0 <- if (ny > 1L) y[nok, , drop = FALSE] else y[nok]
		y  <- if (ny > 1L) y[ ok, , drop = FALSE] else y[ok]
                ## add this information to model.frame as well
                ## need it in outlierStats()
                ## ?? could also add this to na.action, then
                ##    naresid() would pad these as well.
                attr(mf, "zero.weights") <- which(nok)
	    }
	    wts <- sqrt(w)
	    save.y <- y
	    x <- wts * x
	    y <- wts * y
	}
	## check for singular fit

        ## faster, but no longer "allowed" by the Crania:
        ## z0 <- .Call(stats:::C_Cdqrls, x, y, tol = control$solve.tol)
	if(getRversion() >= "3.1.0") {
	    z0 <- .lm.fit(x, y, tol = control$solve.tol)
	    piv <- z0$pivot
	} else {
	    z0 <- lm.fit(x, y, tol = control$solve.tol)
	    piv <- z0$qr$pivot
	}
	rankQR <- z0$rank

	singular.fit <- rankQR < p
	if (rankQR > 0) {
	    if (singular.fit) {
		if (!singular.ok) stop("singular fit encountered")
		pivot <- piv
		p1 <- pivot[seq_len(rankQR)]
		p2 <- pivot[(rankQR+1):p]
		## to avoid problems in the internal fitting methods,
		## split into singular and non-singular matrices,
		## can still re-add singular part later
		dn <- dimnames(x)
		x <- x[,p1]
		attr(x, "assign") <- assign[p1] ## needed for splitFrame to work
	    }
            if (is.function(control$eps.x))
                control$eps.x <- control$eps.x(max(abs(x)))
	    if (!is.null(ini <- init)) {
		if (is.character(init)) {
		    init <- switch(init,
				   "M-S" = lmrob.M.S(x, y, control, mf),
				   "S"   = lmrob.S  (x, y, control, mf=mf),
				   stop('init must be "S", "M-S", function or list'))
		    if(ini == "M-S") { ## "M-S" sometimes reverts to "S":
			ini <- init$control$method
                        ## if(identical(ini, "M-S"))
                        ##     control$method <- paste0(ini, control$method)
                    }
		} else if (is.function(init)) {
		    init <- init(x=x, y=y, control=control, mf=mf)
		} else if (is.list(init)) {
		    ## MK: set init$weights, init$residuals here ??
		    ##	   (needed in lmrob..D..fit)
		    ##	   or disallow method = D... ? would need to fix also
		    ##	  lmrob.kappa: tuning.psi / tuning.chi choice
		    if (singular.fit) {
			## make sure the initial coefficients vector matches
			## to the reduced x
			init$coef <- na.omit(init$coef)
			if (length(init$coef) != ncol(x))
			    stop("Length of initial coefficients vector does not match rank of singular design matrix x")
		    }
		} else stop("unknown init argument")
		stopifnot(is.numeric(init$coef), is.numeric(init$scale))
		## modify (default) control$method, possibly dropping first letter:
		if (control$method == "MM" || substr(control$method, 1, 1) == "S")
		    control$method <- substring(control$method, 2)
		## check for control$cov argument
		if (class(init)[1] != "lmrob.S" && control$cov == '.vcov.avar1')
		    control$cov <- ".vcov.w"
	    }
	    z <- lmrob.fit(x, y, control, init=init, mf = mf) #-> ./lmrob.MM.R
            if(is.character(ini) && !grepl(paste0("^", ini), control$method))
                control$method <- paste0(ini, control$method)
	    if (singular.fit) {
		coef <- numeric(p)
		coef[p2] <- NA
		coef[p1] <- z$coefficients
		names(coef) <- dn[[2L]]
		z$coefficients <- coef
		## Update QR decomposition (z$qr)
		## pad qr and qraux with zeroes (columns that were pivoted to the right in z0)
                d.p <- p-rankQR
                n <- NROW(y)
		z$qr[c("qr","qraux","pivot")] <-
		    list(matrix(c(z$qr$qr, rep.int(0, d.p*n)), n, p,
				dimnames = list(dn[[1L]], dn[[2L]][piv])),
			 ## qraux:
			 c(z$qr$qraux, rep.int(0, d.p)),
			 ## pivot:
			 piv)
	    }
	} else { ## rank 0
	    z <- list(coefficients = if (is.matrix(y)) matrix(NA,p,ncol(y))
				     else rep.int(as.numeric(NA), p),
		      residuals = y, scale = NA, fitted.values = 0 * y,
		      cov = matrix(,0,0), rweights = rep.int(as.numeric(NA), NROW(y)),
		      weights = w, rank = 0, df.residual = NROW(y),
		      converged = TRUE, iter = 0, control=control)
	    if (is.matrix(y)) colnames(z$coefficients) <- colnames(x)
	    else names(z$coefficients) <- colnames(x)
	    if(!is.null(offset)) z$residuals <- y - offset
	}
	if (!is.null(w)) {
	    z$residuals <- z$residuals/wts
	    z$fitted.values <- save.y - z$residuals
	    z$weights <- w
	    if (zero.weights) {
                coef <- z$coefficients
		coef[is.na(coef)] <- 0
		f0 <- x0 %*% coef
		if (ny > 1) {
		    save.r[ok, ] <- z$residuals
		    save.r[nok, ] <- y0 - f0
		    save.f[ok, ] <- z$fitted.values
		    save.f[nok, ] <- f0
		}
		else {
		    save.r[ok] <- z$residuals
		    save.r[nok] <- y0 - f0
		    save.f[ok] <- z$fitted.values
		    save.f[nok] <- f0
		}
		z$residuals <- save.r
		z$fitted.values <- save.f
		z$weights <- save.w
		rw <- z$rweights
		z$rweights <- rep.int(0, length(save.w))
		z$rweights[ok] <- rw
	    }
	}
    }
    if(!is.null(offset))
	z$fitted.values <- z$fitted.values + offset

    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- contrasts
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    z$assign <- assign
    if(control$compute.rd && !is.null(x))
	z$MD <- robMD(x, attr(mt, "intercept"), wqr=z$qr)
    if (model)
	z$model <- mf
    if (ret.x)
	z$x <- if (singular.fit || (!is.null(w) && zero.weights))
	    model.matrix(mt, mf, contrasts) else x
    if (ret.y)
	z$y <- if (!is.null(w)) model.response(mf, "numeric") else y
    class(z) <- "lmrob"
    z
}

if(getRversion() < "3.1.0") globalVariables(".lm.fit")

##' @title Warn about extraneous arguments in the "..."	 (of its caller)
##' @return
##' @author Martin Maechler, June 2012
chk.s <- function(...) {
    if(length(list(...)))
	warning("arguments  ",
		sub(")$", '', sub("^list\\(", '', deparse(list(...), control=c()))),
		"  are disregarded in\n ", deparse(sys.call(-1), control=c()),
		call. = FALSE)
}


##' Robust Mahalanobis Distances
##' internal function, used in lmrob() and plot.lmrob()
robMD <- function(x, intercept, wqr, ...) {
    if(intercept == 1) x <- x[, -1, drop=FALSE]
    if(ncol(x) >= 1) {
	rob <- tryCatch(covMcd(x, ...),
                        warning = function(w) structure("covMcd failed with a warning",
                        class="try-error", condition = w),
                        error = function(e) structure("covMcd failed with an error",
                        class="try-error", condition = e))
	if (inherits(rob, "try-error")) {
            warning("Failed to compute robust Mahalanobis distances, reverting to robust leverages.")
            return(lmrob.leverages(wqr = wqr))
        }
	sqrt( mahalanobis(x, rob$center, rob$cov) )
    }
}

### Method Functions for class lmrob objects ###
### ---------------------------------------- ###

## Many are just wrapper functions for the respective .lm methods

## ---- sorted *ALPHABETICALLY* ----

alias.lmrob <- function(object, ...) {
    ## Purpose: provide alias() for lmrob objects
    ## Cannot use alias.lm directly, since it requires a "clean" object$qr,
    ## i.e., without the robustness weights

    if (is.null(x <- object[["x"]]))
	x <- model.matrix(object)
    weights <- weights(object)
    if (!is.null(weights) && diff(range(weights)))
	x <- x * sqrt(weights)
    object$qr <- qr(x)
    class(object) <- "lm"
    alias(object)
}


## R (3.1.0)-devel copy of case.names.lm() ...../R/src/library/stats/R/lm.R
case.names.lmrob <- function(object, full = FALSE, ...)
{
    w <- weights(object)
    dn <- names(residuals(object))
    if(full || is.null(w)) dn else dn[w!=0]
}

## coef(<lmrob>): no own method ==> using  coef.default(OO) == OO$coefficients
## -------------

## use confint.lm instead of confint.default
## mainly to get t instead of normal quantiles
## Either imported from 'stats' or then copy-paste-defined in ./zzz.R :
confint.lmrob <- confint.lm
dummy.coef.lmrob <- dummy.coef.lm


family.lmrob <- function(object, ...) gaussian() ## == stats:::family.lm


## fitted.default works for "lmrob"

kappa.lmrob <- function(z, ...) kappa.lm(z, ...)

## instead of  stats:::qr.lm()
qrLmr <- function(x) {
    if(!is.list(r <- x$qr))
        stop("lmrob object does not have a proper 'qr' component. Rank zero?")
    r
}

## Basically the same as  stats:::labels.lm -- FIXME: rank 0 fits?
labels.lmrob <- function(object, ...) {
    tl <- attr(object$terms, "term.labels")
    asgn <- object$assign[qrLmr(object)$pivot[seq_len(object$rank)]]
    tl[unique(asgn)]
}

## Works via lm's method [which is still exported]:
model.matrix.lmrob <- model.matrix.lm

## identical to stats:::nobs.lm {but that is hidden .. and small to copy}:
nobs.lmrob <- function(object, ...)
    if (!is.null(w <- object$weights)) sum(w != 0) else NROW(object$residuals)


if(FALSE) ## now replaced with more sophsticated in ./lmrobPredict.R
## learned from MASS::rlm() : via "lm" as well
predict.lmrob <- function (object, newdata = NULL, scale = NULL, ...)
{
    class(object) <- c(class(object), "lm")
    object$qr <- qr(sqrt(object$rweights) * object$x)
    predict.lm(object, newdata = newdata, scale = object$s, ...)
}

print.summary.lmrob <-
    function (x, digits = max(3, getOption("digits") - 3),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"),
              showAlgo = TRUE, ...)
{
    cat("\nCall:\n",
	paste(deparse(x$call, width.cutoff=72), sep = "\n", collapse = "\n"),
	"\n", sep = "")
    control <- lmrob.control.neededOnly(x$control)
    cat(" \\--> method = \"", control$method, '"\n', sep = "")
    ## else cat("\n")
    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    cat(if (!is.null(x$weights) && diff(range(x$weights))) "Weighted ",
	"Residuals:\n", sep = "")
    if (rdf > 5L) {
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <-
	    if (NCOL(resid) > 1)
		structure(apply(t(resid), 1, quantile),
			  dimnames = list(nam, dimnames(resid)[[2]]))
	    else setNames(quantile(resid), nam)
	print(rq, digits = digits, ...)
    }
    else print(resid, digits = digits, ...)
    ## FIXME: need to catch rdf == 0?
    if( length(x$aliased) ) {
	if( !(x$converged) ) {
	    if (x$scale == 0) {
		cat("\nExact fit detected\n\nCoefficients:\n")
	    } else {
		cat("\nAlgorithm did not converge\n")
		if (control$method == "S")
		    cat("\nCoefficients of the *initial* S-estimator:\n")
		else
		    cat(sprintf("\nCoefficients of the %s-estimator:\n",
				control$method))
	    }
	    printCoefmat(x$coef, digits = digits, signif.stars = signif.stars,
			 ...)
	} else {
	    if (nsingular <- df[3L] - df[1L])
		cat("\nCoefficients: (", nsingular,
		    " not defined because of singularities)\n", sep = "")
	    else cat("\nCoefficients:\n")
	    coefs <- x$coefficients
	    if(!is.null(aliased <- x$aliased) && any(aliased)) {
		cn <- names(aliased)
		coefs <- matrix(NA, length(aliased), 4, dimnames=list(cn, colnames(coefs)))
		coefs[!aliased, ] <- x$coefficients
	    }

	    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
			 na.print="NA", ...)
	    cat("\nRobust residual standard error:",
		format(signif(x$scale, digits)),"\n")
          if (!is.null(x$r.squared) && x$df[1] != attr(x$terms, "intercept")) {
                cat("Multiple R-squared: ", formatC(x$r.squared, digits = digits))
                cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared, digits = digits),
                    "\n")
            }
	    ## FIXME: use naprint() here to list observations deleted due to missingness?
	    correl <- x$correlation
	    if (!is.null(correl)) {
		p <- NCOL(correl)
		if (p > 1) {
		    cat("\nCorrelation of Coefficients:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			print(symnum(correl), abbr.colnames = NULL)
		    }
		    else { correl <- format(round(correl, 2), nsmall = 2,
					    digits = digits)
			   correl[!lower.tri(correl)] <- ""
			   print(correl[-1, -p, drop = FALSE], quote = FALSE)
		       }
		}
	    }
	    cat("Convergence in", x$iter, "IRWLS iterations\n")
	}
	cat("\n")

	if (!is.null(rw <- x$rweights)) {
	    if (any(zero.w <- x$weights == 0))
		rw <- rw[!zero.w]
            eps.outlier <- if (is.function(EO <- control$eps.outlier))
                EO(nobs(x)) else EO
	    summarizeRobWeights(rw, digits = digits, eps = eps.outlier, ...)
	}

    } else cat("\nNo Coefficients\n")

    if (showAlgo && !is.null(control))
	printControl(control, digits = digits, drop. = "method")
    invisible(x)
}


print.lmrob <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n", cl <- deparse(x$call, width.cutoff=72), "\n", sep = "")
    control <- lmrob.control.neededOnly(x$control)
    if(!any(grepl("method *= *['\"]", cl)))## 'method = ".."' not explicitly visible above
	cat(" \\--> method = \"", control$method, '"\n', sep = "") else cat("\n")
    if(length((cf <- coef(x)))) {
	if( x$converged )
	    cat("Coefficients:\n")
	else {
	    if (x$scale == 0) {
		cat("Exact fit detected\n\nCoefficients:\n")
	    } else {
		cat("Algorithm did not converge\n\n")
		if (control$method == "S")
		    cat("Coefficients of the *initial* S-estimator:\n")
		else
		    cat(sprintf("Coefficients of the %s-estimator:\n",
				control$method))
	    }
	}
	print(format(cf, digits = digits), print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

print.lmrob.S <- function(x, digits = max(3, getOption("digits") - 3),
			  showAlgo = TRUE, ...)
{
    cat("S-estimator lmrob.S():\n")
    if(length((cf <- coef(x)))) {
	if (x$converged)
	    cat("Coefficients:\n")
	else if (x$scale == 0)
	    cat("Exact fit detected\n\nCoefficients:\n")
	else
	    cat("Algorithm did not converge\n\n")
	print(format(cf, digits = digits), print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n")
    cat("scale = ",format(x$scale, digits=digits), "; ",
	if(x$converged)"converged" else "did NOT converge",
	" in ", x$k.iter, " refinement steps\n")
    if (showAlgo && !is.null(x$control))
        printControl(x$control, digits = digits, drop. = "method")
    invisible(x)
}


## practically identical to  stats:::qr.lm :
qr.lmrob <- function (x, ...) {
    if (is.null(r <- x$qr))
	stop("lmrob object does not have a proper 'qr' component. Rank must be zero")
    r
}

residuals.lmrob <- function(object, ...) residuals.lm(object, ...)

## even simpler than residuals.default():
residuals.lmrob.S <- function(obj) obj$residuals

summary.lmrob <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    if (is.null(object$terms))
	stop("invalid 'lmrob' object:  no terms component")
    p <- object$rank
    df <- object$df.residual #was $degree.freedom
    sigma <- object[["scale"]]
    aliased <- is.na(coef(object))
    cf.nms <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    if (p > 0) {
	n <- p + df
	p1 <- seq_len(p)
	se <- sqrt(if(length(object$cov) == 1L) object$cov else diag(object$cov))
	est <- object$coefficients[object$qr$pivot[p1]]
	tval <- est/se
	ans <- object[c("call", "terms", "residuals", "scale", "rweights",
			"converged", "iter", "control")]
	if (!is.null(ans$weights))
	    ans$residuals <- ans$residuals * sqrt(object$weights)
	## 'df' vector, modeled after summary.lm() : ans$df <- c(p, rdf, NCOL(Qr$qr))
	## where  p <- z$rank ; rdf <- z$df.residual ; Qr <- qr.lm(object)
	ans$df <- c(p, df, NCOL(object$qr$qr))
	ans$coefficients <-
	    if( ans$converged)
		cbind(est, se, tval, 2 * pt(abs(tval), df, lower.tail = FALSE))
	    else
		cbind(est, if(sigma <= 0) 0 else NA, NA, NA)
	dimnames(ans$coefficients) <- list(names(est), cf.nms)
        if (p != attr(ans$terms, "intercept")) {
            df.int <- if (attr(ans$terms, "intercept")) 1L else 0L
            ## This block is based on code by Olivier Renaud <Olivier.Renaud@unige.ch>
            resid <- object$residuals
            pred <- object$fitted.values
            resp <- if (is.null(object[["y"]])) pred + resid else object$y
            wgt <- object$rweights
            ## scale.rob <- object$scale
            ## correction E[wgt(r)] / E[psi'(r)] ( = E[wgt(r)] / E[r*psi(r)] )
            ctrl <- object$control
            c.psi <- ctrl$tuning.psi
            psi <- ctrl$psi
            correc <-
                if (psi == 'ggw') {
                    if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.95, NA)))) 1.121708
                    else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.95, NA)))) 1.163192
                    else if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.85, NA)))) 1.33517
                    else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.85, NA)))) 1.395828
                    else lmrob.E(wgt(r), ctrl) / lmrob.E(r*psi(r), ctrl)
		} else if (any(psi == .Mpsi.R.names) &&
			   isTRUE(all.equal(c.psi, .Mpsi.tuning.default(psi)))) {
                    switch(psi,
                           bisquare = 1.207617,
                           welsh    = 1.224617,
                           optimal  = 1.068939,
                           hampel   = 1.166891,
                           lqq      = 1.159232,
			   stop('unsupported psi function -- should not happen'))
                } else lmrob.E(wgt(r), ctrl) / lmrob.E(r*psi(r), ctrl)
            resp.mean <- if (df.int == 1L) sum(wgt * resp)/sum(wgt) else 0
            yMy <- sum(wgt * (resp - resp.mean)^2)
            rMr <- sum(wgt * resid^2)
            ans$r.squared <- r2correc <- (yMy - rMr) / (yMy + rMr * (correc - 1))
            ans$adj.r.squared <- 1 - (1 - r2correc) * ((n - df.int) / df)
        } else ans$r.squared <- ans$adj.r.squared <- 0
	ans$cov <- object$cov
	if(length(object$cov) > 1L)
	    dimnames(ans$cov) <- dimnames(ans$coefficients)[c(1,1)]
	if (correlation) {
	    ans$correlation <- ans$cov / outer(se, se)
	    ans$symbolic.cor <- symbolic.cor
	}
    } else { ## p = 0: "null model"
	ans <- object
	ans$df <- c(0L, df, length(aliased))
	ans$coefficients <- matrix(NA, 0L, 4L, dimnames = list(NULL, cf.nms))
        ans$r.squared <- ans$adj.r.squared <- 0
	ans$cov <- object$cov
    }
    ans$aliased <- aliased # used in print method
    ans$sigma <- sigma # 'sigma': in summary.lm() & 'fit.models' pkg
    if (is.function(ans$control$eps.outlier))
        ans$control$eps.outlier <- ans$control$eps.outlier(nobs(object))
    if (is.function(ans$control$eps.x))
	ans$control$eps.x <-
	    if(!is.null(o.x <- object[['x']]))
                ans$control$eps.x(max(abs(o.x))) ## else NULL
    structure(ans,
	      class = "summary.lmrob")
}


## R (3.1.0)-devel copy of variable.names.lm() ...../R/src/library/stats/R/lm.R
variable.names.lmrob <- function(object, full = FALSE, ...)
{
    if(full) dimnames(qrLmr(object)$qr)[[2L]]
    else if(object$rank) dimnames(qrLmr(object)$qr)[[2L]][seq_len(object$rank)]
    else character()
}

vcov.lmrob <- function (object, cov = object$control$cov, ...) {
    if (!is.null(object$cov) && (missing(cov) ||
				 identical(cov, object$control$cov)))
	object$cov
    else {
	## cov = ..$control$cov is typically ".vcov.w" or ".vcov.avar1"
	lf.cov <- if (!is.function(cov)) get(cov, mode = "function") else cov
	lf.cov(object, ...)
    }
}

sigma.lmrob <- function(object, ...) object$scale

weights.lmrob <- function(object, type = c("prior", "robustness"), ...) {
    type <- match.arg(type)
    res <- if (type == "prior") {
	## Issue warning only if called from toplevel. Otherwise the warning pop
	## up at quite unexpected places, e.g., case.names().
	if (is.null(object[["weights"]]) && identical(parent.frame(), .GlobalEnv))
	    warning("No weights defined for this object. Use type=\"robustness\" argument to get robustness weights.")
	object[["weights"]]
    } else object[["rweights"]]
    if (is.null(object$na.action))
	res
    else naresid(object$na.action, res)
}


####  functions hidden in namespace ####

printControl <-
    function(ctrl, digits = getOption("digits"),
	     str.names = "seed", drop. = character(0),
	     header = "Algorithmic parameters:",
	     ...)
{
    ## Purpose: nicely and sensibly print a 'control' structure
    ##		currently  for lmrob(), glmrob()
    ## Author: Martin Maechler, Date: 31 May 2006
    PR <- function(LST, ...) if(length(LST)) print(unlist(LST), ...)

    cat(header,"\n")
    is.str <- (nc <- names(ctrl)) %in% str.names
    do. <- !is.str & !(nc %in% drop.)
    is.ch <- vapply(ctrl, is.character, NA)
    real.ctrl <- vapply(ctrl, function(x)
			length(x) > 0 && is.numeric(x) && any(x %% 1 != 0), NA)
    PR(ctrl[do. & real.ctrl], digits = digits, ...)
    ## non-real, non-char ones (typically integers), but dropping 0-length ones
    PR(ctrl[do. & !is.ch & !real.ctrl], ...)
    ## char ones
    PR(ctrl[do. & is.ch], ...)
    if(any(is.str))
	for(n in nc[is.str]) {
	    cat(n,":")
	    str(ctrl[[n]], vec.len = 2)
	    ## 'vec.len = 2' is smaller than normal, but nice for Mersenne seed
	}
}


summarizeRobWeights <-
    function(w, digits = getOption("digits"), header = "Robustness weights:",
	     eps = 0.1 / length(w), eps1 = 1e-3, ...)
{
    ## Purpose: nicely print a "summary" of robustness weights
    stopifnot(is.numeric(w))
    cat(header,"\n")
    cat0 <- function(...) cat('', ...)
    n <- length(w)
    if(n <= 10) print(w, digits = digits, ...)
    else {
	n1 <- sum(w1 <- abs(w - 1) < eps1)
	n0 <- sum(w0 <- abs(w) < eps)
	if(any(w0 & w1))
	    warning("weights should not be both close to 0 and close to 1!\n",
		    "You should use different 'eps' and/or 'eps1'")
	if(n0 > 0 || n1 > 0) {
	    if(n0 > 0) {
		formE <- function(e) formatC(e, digits = max(2, digits-3), width=1)
		i0 <- which(w0)
		maxw <- max(w[w0])
		c3 <- paste0("with |weight| ",
                             if(maxw == 0) "= 0" else paste("<=", formE(maxw)),
			    " ( < ", formE(eps), ");")
		cat0(if(n0 > 1) {
		       cc <- sprintf("%d observations c(%s)",
				     n0, strwrap(paste(i0, collapse=",")))
		       c2 <- " are outliers"
		       paste0(cc,
			     if(nchar(cc)+ nchar(c2)+ nchar(c3) > getOption("width"))
			     "\n	", c2)
		     } else
		       sprintf("observation %d is an outlier", i0),
		     c3, "\n")
	    }
	    if(n1 > 0)
		cat0(ngettext(n1, "one weight is",
			     sprintf("%s%d weights are",
				     if(n1 == n)"All " else '', n1)), "~= 1.")
	    n.rem <- n - n0 - n1
	    if(n.rem <= 0) { # < 0 possible if w0 & w1 overlap
		if(n1 > 0) cat("\n")
		return(invisible())
	    }
	    cat0("The remaining",
		 ngettext(n.rem, "one", sprintf("%d ones", n.rem)), "are")
	    if(is.null(names(w)))
		names(w) <- as.character(seq(along = w))
	    w <- w[!w1 & !w0]
	    if(n.rem <= 10) {
		cat("\n")
		print(w, digits = digits, ...)
		return(invisible())
	    }
	    else cat(" summarized as\n")
	}
	print(summary(w, digits = digits), digits = digits, ...)
    }
}
