#### This is originally from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov

##  I would like to thank Peter Rousseeuw and Katrien van Driessen for
##  providing the initial code of this function.

### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, a copy is available at
### http://www.r-project.org/Licenses/


ltsReg <- function(x, ...) UseMethod("ltsReg")

ltsReg.formula <- function(formula, data, subset, weights, na.action,
			   model = TRUE, x.ret = FALSE, y.ret = FALSE,
                           contrasts = NULL, offset, ...)
{
    cl <- match.call()
    ##	  method <- match.arg(method)

    ## keep only the arguments which should go into the model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    ##	  if (method == "model.frame") return(mf)

    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric") ## was model.extract(mf, "response")

    if (is.empty.model(mt)) { # "y ~ 0" : no coefficients
	x <- offset <- NULL
	fit <- list(method = "ltsReg for empty model",
		    coefficients = numeric(0), residuals = y,
		    fitted.values = 0 * y, lts.wt = 1 + 0 * y,
		    rank = 0, intercept = FALSE, df.residual = length(y))
	## alpha = alpha from "..."
	class(fit) <- "lts"
    }
    else {
        w <- model.weights(mf)
        offset <- model.offset(mf)

	x <- model.matrix(mt, mf, contrasts)

	## Check if there is an intercept in the model.
	## A formula without intercept looks like this: Y ~ . -1
	## If so, remove the corresponding column and use intercept=TRUE
	## in the call to ltsReg.default(); by default, intercept=FALSE.
	xint <- match("(Intercept)", colnames(x), nomatch = 0)
	if(xint)
	    x <- x[, -xint, drop = FALSE]
	fit <- ltsReg.default(x, y, intercept = (xint > 0), ...)
    }

    ## 3) return the na.action info
    fit$na.action <- attr(mf, "na.action")
    fit$offset <- offset

    ## 4) return the contrasts used in fitting: possibly as saved earlier.
    fit$contrasts <- attr(x, "contrasts")

    fit$xlevels <- .getXlevels(mt, mf)
    fit$call <- cl
    fit$terms <- mt

    if(model) fit$model <- mf
    if(x.ret) fit$x <- x # or? if(xint == 0) x else  x[, c(2:p,1), drop=FALSE]
    if(y.ret) fit$y <- y

    fit
}

ltsReg.default <- function (x, y, intercept = TRUE,
        alpha = control$ alpha,
        nsamp = control$ nsamp,
	adjust = control$ adjust,
	mcd = TRUE,
	qr.out = FALSE,
	yname = NULL,
        seed  = control$ seed,
        trace = control$ trace,
	use.correction = control$ use.correction,
        wgtFUN = control$ wgtFUN,
	control = rrcov.control(),
	...)
{
    ##	 Analyze and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.
### MM: FIXME: this sucks ('control' may contain *some* but not all parts!):
    if(!missing(control)) {
	defCtrl <- rrcov.control()	# default control
	if(is.null(alpha) && control$alpha != defCtrl$alpha)
	    alpha <- control$alpha
	if(nsamp == defCtrl$nsamp)    nsamp <- control$nsamp
	if(identical(seed, defCtrl$seed)) seed <- control$seed

	if(use.correction == defCtrl$use.correction)
	    use.correction <- control$use.correction
	if(adjust == defCtrl$adjust)
	    adjust <- control$adjust
    } else defCtrl <- control ## == rrcov.control()
    ## For back compatibility, as some new args did not exist pre 2013-04,
    ## and callers of covMcd() may use a "too small"  'control' list:

    if(missing(wgtFUN)) getDefCtrl("wgtFUN", defCtrl)

    if(length(seed) > 0) {
	if(length(seed) < 3 || seed[1L] < 100)
	    stop("invalid 'seed'. Must be compatible with .Random.seed !")
	if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
	    seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
	    on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
	}
	assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    if(alpha < 1/2 || alpha > 1)
	stop("alpha not inside [1/2, 1]")

    ## FIXME: change this analogously to covMcd()'s and covComedian()'s
    ## quantiel <- qnorm(0.9875)
    if(is.character(wgtFUN)) {
	switch(wgtFUN,
	       "01.original" = {
		   cW <- qnorm(0.9875)
		   wgtFUN <- function(r) as.numeric(abs(r) <= cW)
	       },
	       stop("unknown 'wgtFUN' specification: ", wgtFUN))
    } else if(!is.function(wgtFUN))
	stop("'wgtFUN' must be a function or a string specifying one")

    ## vt::03.02.2006 - raw.cnp2 and cnp2 are vectors of size 2 and  will
    ##	 contain the correction factors (concistency and finite sample)
    ##	 for the raw and reweighted estimates respectively. Set them initially to 1.
    ##	 If use.correction is set to FALSE (default=TRUE), the finite sample correction
    ##	 factor will not be used (neither for the raw estimates nor for the reweighted)
    raw.cnp2 <- rep(1,2)
    cnp2 <- rep(1,2)

    ##cat("++++++ Entering ltsReg() ...\n")

    y <- data.matrix(y)
    if (!is.numeric(y)) stop("y is not a numeric")
    if (dim(y)[2] != 1) stop("y is not onedimensional")

    oneD <- (missing(x) || is.null(x) || NCOL(x) == 0) ## location model - no x
    if(oneD) {
	x <- matrix(1, nrow(y), 1)
    }
    else { ## x is present
	if(is.data.frame(x))
	    x <- data.matrix(x)
	else if (!is.matrix(x))
	    x <- matrix(x, length(x), 1,
			dimnames = list(names(x), deparse(substitute(x))))
    }

    if (nrow(x) != nrow(y))
	stop("Number of observations in x and y not equal")
    na.x <- !is.finite(rowSums(x))
    na.y <- !is.finite(y)
    ok <- !(na.x | na.y)
    x <- x[ok, , drop = FALSE]
    y <- y[ok, , drop = FALSE]
    dx <- dim(x)
    n <- dx[1]
    if (n == 0)
	stop("All observations have missing values!")
    dimny <- dimnames(y)
    rownames <- dimny[[1]]
    yn <- if(!is.null(yname))
	yname else if(!is.null(dimny[[2]])) dimny[[2]]
    has.yn <- !is.null(yn)
    if(!has.yn) yn <- "Y"
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    if (!oneD) {
	is.const <- function(x) {
	    c1 <- range(x)
	    c1[1] == c1[2]
	}
	if (any(apply(x, 2, is.const)))
	    stop("There is at least one constant column. Remove it and set intercept=TRUE")
    }

    ##cat("++++++ Prepare: Ready.\n")

    xn <- (dnx <- dimnames(x))[[2]]
    xn <- if(!is.null(xn)) xn else if (dx[2] > 1)
	paste("X", 1:dx[2], sep = "") else if (dx[2]) "X" ## else : p = 0
    dimnames(x) <- list(dnx[[1]], xn) # also works  if(is.null(dnx))
    y <- as.vector(y)

    if(all(x == 1)) { ## includes 'oneD' and empty x (p = 0)
	if(qr.out) {
	    warning("'qr.out = TRUE' for univariate location is disregarded")
	    qr.out <- FALSE
	}
	h <- h.alpha.n(alpha, n, dx[2])
	p <- 1
	if (alpha == 1) {
	    scale <- sqrt(cov.wt(as.matrix(y))$cov)
	    center <- as.vector(mean(y))
	    ## xbest <- NULL
	} else {
            sh <- .fastmcd(as.matrix(y), as.integer(h), nsamp = 0, # (y *is* 1-dim.!)
			   nmini = 300, kmini = 5)
	    center <- as.double(sh$initmean)
	    qalpha <- qchisq(h/n, 1)
	    calphainvers <- pgamma(qalpha/2, 1/2 + 1)/(h/n)
	    raw.cnp2[1] <- calpha <- 1/calphainvers
	    raw.cnp2[2] <- correct <- LTScnp2(1, intercept = intercept, n, alpha)
	    if(!use.correction) # do not use finite sample correction factor
		raw.cnp2[2] <- correct <- 1.0

	    scale <- sqrt(as.double(sh$initcovariance)) * sqrt(calpha) * correct
	    ## xbest <- sort(as.vector(sh$best))  # fastmcd in the univariate case does not return inbest[]
	}
	resid <- y - center
	ans <- list(method = "Univariate location and scale estimation.",
		    best = NULL, # xbest,
		    coefficients = center,
		    alpha = alpha,
		    quan  = h,
		    raw.coefficients = center,
		    raw.resid = resid/scale,
		    raw.weights = rep.int(NA, length(na.y)))
	if(abs(scale) < 1e-07) {
	    ans$raw.weights[ok] <- weights <- as.numeric(abs(resid) < 1e-07)
	    ans$scale <- ans$raw.scale <- 0
	    ans$crit <- 0
	    ans$method <- paste(ans$method,
				"More than half of the data are equal!",sep="\n")
	}
	else {
	    ans$raw.scale <- scale
	    ans$raw.weights[ok] <- weights <- wgtFUN(resid/scale)
            sum.w <- sum(weights)
	    reweighting <- cov.wt(as.matrix(y), wt = weights)
	    ans$coefficients <- reweighting$center
	    ans$scale <- sqrt(sum.w/(sum.w - 1) * reweighting$cov)
	    resid <- y - ans$coefficients
	    ans$crit <- sum(sort((y - center)^2, partial = h)[1:h])
	    if (sum.w != n) {
		qdelta.rew <- qchisq(sum.w/n, 1)
		cdeltainvers.rew <- pgamma(qdelta.rew/2, 1/2 + 1)/(sum.w/n)
		cdelta.rew <- sqrt(1/cdeltainvers.rew)
		correct.rew <-
		    if(use.correction)
			LTScnp2.rew(1, intercept = intercept, n, alpha) else 1
		cnp2 <- c(cdelta.rew, correct.rew)
		ans$scale <- ans$scale * cdelta.rew * correct.rew
	    }
	    weights <- wgtFUN(resid/ans$scale)
	}
	fitted <- ans$coefficients
	ans$resid <- resid/ans$scale
	ans$rsquared <- 0
	ans$intercept <- intercept
        if(has.yn)
            names(ans$coefficients) <- names(ans$raw.coefficients) <- yn

    } ## end {all(x == 1)} --

    else { ## ------------------ usual non-trivial case ---------------------
	if(mcd) ## need 'old x' later
	    X <- x
	if (intercept) { ## intercept must be *last* (<- fortran code) {"uahh!"}
	    x <- cbind(x, "Intercept" = 1)
	    dx <- dim(x)
	    xn <- colnames(x)
	}
	p <- dx[2]
	if (n <= 2 * p)
	    stop("Need more than twice as many observations as variables.")

	## VT:: 26.12.2004
	## Reorder the coefficients so that the intercept is at the beginning ..
	getCoef <- ## simple wrapper (because of above "intercept must be")
	    if(p > 1 && intercept)
		 function(cf) cf[c(p, 1:(p - 1))]
	    else function(cf) cf

	ans <- list(alpha = alpha, raw.weights = rep.int(NA, length(na.y)))

	if(alpha == 1) { ## alpha == 1 -----------------------
	    ## old, suboptimal: z <- lsfit(x, y, intercept = FALSE)
	    z <- lm.fit(x, y)
	    qrx <- z$qr
	    cf <- z$coef
	    names(cf) <- xn
	    ans$raw.coefficients <- getCoef(cf)

	    resid <- z$residuals
	    ans$quan <- h <- n

	    s0 <- sqrt((1/(n - p)) * sum(resid^2))

	    ##cat("++++++ B - alpha == 1... - s0=",s0,"\n")
	    if(abs(s0) < 1e-07) {
		fitted <- x %*% z$coef
		ans$raw.weights[ok] <- weights <- as.numeric(abs(resid) <= 1e-07)
		ans$scale <- ans$raw.scale <- 0
		ans$coefficients <- ans$raw.coefficients
	    }
	    else {
		ans$raw.scale <- s0
		ans$raw.resid <- resid / s0
		ans$raw.weights[ok] <- weights <- wgtFUN(ans$raw.resid)
		sum.w <- sum(weights)
		## old, suboptimal: z <- lsfit(x, y, wt = weights, intercept = FALSE)
		z <- lm.wfit(x, y, w = weights)

		ans$coefficients <- getCoef(z$coef)

		fitted <- x %*% z$coef
		ans$scale <- sqrt(sum(weights * resid^2)/(sum.w - 1))
		if (sum.w != n) {
		    qn.w <- qnorm((sum.w + n)/(2 * n))
		    cdelta.rew <- 1/sqrt(1 - (2 * n)/(sum.w/qn.w) * dnorm(qn.w))
		    ans$scale <- ans$scale * cdelta.rew
		}
		ans$resid <- resid/ans$scale
		weights <- wgtFUN(ans$resid)
	    }

	    names(ans$coefficients) <- getCoef(xn)

	    s1 <- sum(resid^2)
	    ans$crit <- s1
	    sh <- (if (intercept) y - mean(y) else y) ^ 2
	    ans$rsquared <- max(0, min(1, 1 - (s1/sh)))

	    ans$method <- "Least Squares Regression."

	} ## end {alpha == 1} : "classical"

	else { ## alpha < 1 -----------------------------------------------
	    coefs <- rep(NA, p)
	    names(coefs) <- xn
	    qrx <- if(qr.out) qr(x) else qr(x)[c("rank", "pivot")]

	    rk <- qrx$rank
	    if (rk < p)
		stop("x is singular")
	    ## else :

	    h <- h.alpha.n(alpha, n, rk)

	    z <- .fastlts(x, y, h, nsamp, intercept, adjust, trace=as.integer(trace))
	    if(z$objfct < 0)
		stop("no valid subsample found in LTS - set 'nsamp' or rather use lmrob.S()")
	    ## vt:: lm.fit.qr == lm.fit(...,method=qr,...)
	    cf <- lm.fit(x[z$inbest, , drop = FALSE], y[z$inbest])$coef
	    if(any(ic <- is.na(cf)))
		stop(gettextf("NA coefficient (at %s) from \"best\" subset",
			      paste(which(ic), collapse =",")))
	    ans$best <- sort(z$inbest)
	    fitted <- x %*% cf
	    resid <- y - fitted
	    piv <- 1:p
	    coefs[piv] <- cf ## FIXME? why construct 'coefs' so complicatedly?	use 'cf' !

	    ans$raw.coefficients <- getCoef(coefs)

	    ans$quan <- h
	    correct <- if(use.correction)
		LTScnp2(p, intercept = intercept, n, alpha) else 1
	    raw.cnp2[2] <- correct
	    s0 <- sqrt(mean(sort(resid^2, partial = h)[1:h]))
	    sh0 <- s0
	    qn.q <- qnorm((h + n)/ (2 * n))
	    s0 <- s0 / sqrt(1 - (2 * n)/(h / qn.q) * dnorm(qn.q)) * correct

	    if (abs(s0) < 1e-07) {
		ans$raw.weights[ok] <- weights <- as.numeric(abs(resid) <= 1e-07)
		ans$scale <- ans$raw.scale <- 0
		ans$coefficients <- ans$raw.coefficients
	    }
	    else {
		ans$raw.scale <- s0
		ans$raw.resid <- resid/ans$raw.scale
		ans$raw.weights[ok] <- weights <- wgtFUN(resid/s0)
		sum.w <- sum(weights)

		## old, suboptimal: z1 <- lsfit(x, y, wt = weights, intercept = FALSE)
		z1 <- lm.wfit(x, y, w = weights)

		ans$coefficients <- getCoef(z1$coef)

		fitted <- x %*% z1$coef
		resid <- z1$residuals
		ans$scale <- sqrt(sum(weights * resid^2)/(sum.w - 1))
		if (sum.w == n) {
		    cdelta.rew <- 1
		    correct.rew <- 1
		}
		else {
		    qn.w <- qnorm((sum.w + n)/(2 * n))
		    cnp2[1] <- cdelta.rew <- 1 / sqrt(1 - (2 * n)/(sum.w / qn.w) * dnorm(qn.w))
		    correct.rew <-
			if (use.correction) ## use finite sample correction
			    LTScnp2.rew(p, intercept = intercept, n, alpha)
			else 1
		    cnp2[2] <- correct.rew
		    ans$scale <- ans$scale * cdelta.rew * correct.rew
		}
		ans$resid <- resid/ans$scale
		weights <- wgtFUN(ans$resid)
	    }
	    ## unneeded: names(ans$coefficients) <- names(ans$raw.coefficients)
	    ans$crit <- z$objfct
	    if (intercept) {
                sh <- .fastmcd(as.matrix(y), as.integer(h), nsamp = 0, # (y *is* 1-dim.!)
			       nmini = 300, kmini = 5)
		y <- as.vector(y) ## < ??
		sh <- as.double(sh$adjustcov)
		iR2 <- (sh0/sh)^2
	    }
	    else {
		s1 <- sum(sort(resid^2, partial = h)[1:h])
		sh <- sum(sort(y^2,     partial = h)[1:h])
		iR2 <- s1/sh
	    }

	    ans$rsquared <- if(is.finite(iR2)) max(0, min(1, 1 - iR2)) else 0

	    attributes(resid) <- attributes(fitted) <- attributes(y)
	    ans$method <- "Least Trimmed Squares Robust Regression."
	} ## end { alpha < 1 }

	ans$intercept <- intercept
	if (abs(s0) < 1e-07)
	    ans$method <- paste(ans$method, "\nAn exact fit was found!")

	if (mcd) { ## compute robust distances {for diagnostics, eg. rdiag()plot}
	    mcd <- covMcd(X, alpha = alpha, use.correction=use.correction)
	    if ( -determinant(mcd$cov, logarithm = TRUE)$modulus > 50 * p) {
		ans$RD <- "singularity"
	    }
	    else {
		ans$RD <- rep.int(NA, length(na.y))
		ans$RD[ok] <- sqrt(mahalanobis(X, mcd$center, mcd$cov))
		names(ans$RD) <- rownames
	    }
	}

    } ## end { nontrivial 'x' }

    ans$lts.wt <- rep.int(NA, length(na.y))
    ans$lts.wt[ok] <- weights
    ans$residuals <- rep.int(NA, length(na.y))
    ans$residuals[ok] <- resid
    ans$fitted.values <- rep.int(NA, length(na.y))
    ans$fitted.values[ok] <- fitted

    names(ans$fitted.values) <- names(ans$residuals) <- names(ans$lts.wt) <-
	rownames
    if(has.yn) { ## non-sense otherwise:
	names(ans$scale) <- names(ans$raw.scale) <- yn
	names(ans$rsquared) <- names(ans$crit) <- yn
    }
    ans$Y <- y
    ans$X <- if(p > 1 && intercept) x[, c(p, 1:(p - 1))] else x
    dimnames(ans$X) <- list(rownames[ok], names(ans$coefficients))
    if (qr.out)
	ans$qr <- qrx
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    class(ans) <- "lts"
    ans$call <- match.call()
    return(ans)
} ## {ltsReg.default}

summary.lts <- function (object, correlation = FALSE, ...)
{
    z <- object
    r <- z$residuals
    f <- z$fitted
    int <- z$intercept
    w <- as.vector(z$lts.wt)
    n <- sum(w)

    Qr <- qr(w * z$X)# 'w * z$X': more efficient than t(t(object$X) %*% diag(w))
    p <- Qr$rank
    p1 <- seq(length = p) ## even for p = 0
    rdf <- n - p
    mss <-  if(int) {
		m <- sum(w * f /sum(w))
		sum(w * (f - m)^2)
	    } else
		sum(w * f^2)
    rss <- sum(w * r^2)

    r <- sqrt(w) * r
    resvar <- rss/rdf

    R <- if (p > 0) chol2inv(Qr$qr[p1, p1, drop = FALSE]) else matrix(,p,p)
    ## no need to reorder R anymore, since 'X' already has "intercept first"
    se <- sqrt(diag(R) * resvar)

    est <- z$coefficients
    tval <- est/se

    ans <-
	c(z[c("call", "terms")],
	  ## not again attr(ans, "call") <- attr(z,"call")
	  list(residuals = r,
	       coefficients = {
		   cbind("Estimate" = est, "Std. Error" = se, "t value" = tval,
			 "Pr(>|t|)" = 2*pt(abs(tval), rdf, lower.tail = FALSE))
	       },
	       sigma = sqrt(resvar),
	       df = c(p, rdf, NCOL(Qr$qr))))

    df.int <- if(int) 1 else 0
    if(p - df.int > 0) {
	ans$r.squared <- mss/(mss + rss)
	ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
	ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
			    numdf = p - df.int, dendf = rdf)
    } else
	ans$r.squared <- ans$adj.r.squared <- 0

    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]

    if (correlation) {
	ans$correlation <- (R * resvar)/outer(se, se)
	dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    }
    class(ans) <- "summary.lts"
    ans
}

print.lts <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
	cat("Coefficients:\n")
	print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
	cat("\nScale estimate", format(x$scale, digits = digits) ,"\n\n")
    }
    else
	cat("No coefficients\n")

    invisible(x)
}

print.summary.lts <-
    function(x, digits = max(3, getOption("digits") - 3),
	     signif.stars = getOption("show.signif.stars"), ...)
##	     signif.stars = FALSE, ...)
##			    ^^^^^ (since they are not quite correct ?)
{
    cat("\nCall:\n",
	paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    resid <- x$residuals
    df <- x$df
    rdf <- df[2]
    cat("Residuals (from reweighted LS):\n")
    ## "cut & paste" from print.summary.lm():
    if(rdf > 5) {
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <-	if(length(dim(resid)) == 2)
		    structure(apply(t(resid), 1, quantile),
			      dimnames = list(nam, dimnames(resid)[[2]]))
		else
		    structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
    }
    else if(rdf > 0) {
	print(resid, digits = digits, ...)
    } else { # rdf == 0 : perfect fit!
	cat("ALL", df[1], "residuals are 0: no residual degrees of freedom!\n")
    }

    if(NROW(x$coefficients)) {
	if (nsingular <- df[3] - df[1])
	    cat("\nCoefficients: (", nsingular,
		" not defined because of singularities)\n", sep = "")
	else
	    cat("\nCoefficients:\n")
	printCoefmat(x$coefficients, digits = digits,
		     signif.stars = signif.stars, ...)
    }
    else cat("\nNo coefficients\n")

    cat("\nResidual standard error:",
    format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")

    if(!is.null(x$fstatistic)) {
	cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
	cat(",\tAdjusted R-squared:",formatC(x$adj.r.squared,digits = digits),
	    "\nF-statistic:", formatC(x$fstatistic[1], digits = digits),
	    "on", x$fstatistic[2], "and",
	    x$fstatistic[3], "DF,  p-value:",
	    format.pval(pf(x$fstatistic[1], x$fstatistic[2],
			   x$fstatistic[3], lower.tail = FALSE), digits = digits),
	    "\n")
    }

    correl <- x$correlation
    if(!is.null(correl)) {
	p <- NCOL(correl)
	if(p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    correl <- format(round(correl, 2), nsmall = 2, digits = digits)
	    correl[!lower.tri(correl)] <- ""
	    print(correl[-1, -p, drop = FALSE], quote = FALSE)
	}
    }
    cat("\n")
    invisible(x)
}


### --- Namespace hidden (but parsed once and for all) : -------------

##' Compute Finite Sample Correction Factor for the "raw" LTSreg() scale
LTScnp2 <- function(p, intercept = intercept, n, alpha)
{
    stopifnot(0.5 <= alpha, alpha <= 1)
    if (intercept)
	p <- p - 1
    stopifnot(p == as.integer(p), p >= 0)
    if (p == 0) {
	fp.500.n <- 1 - exp( 0.262024211897096) / n^ 0.604756680630497
	fp.875.n <- 1 - exp(-0.351584646688712) / n^ 1.01646567502486
	if ((0.5 <= alpha) && (alpha <= 0.875)) {
	    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
	    fp.alpha.n <- sqrt(fp.alpha.n)
	}
	if ((0.875 < alpha) && (alpha < 1)) {
	    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)
	    fp.alpha.n <- sqrt(fp.alpha.n)
	}
    }
    else { ## p >= 1
	if (p == 1) {
	    if (intercept) {
		fp.500.n <- 1 - exp( 0.630869217886906 ) / n^ 0.650789250442946
		fp.875.n <- 1 - exp( 0.565065391014791 ) / n^ 1.03044199012509
	    }
	    else {
		fp.500.n <- 1 - exp(-0.0181777452315321) / n^ 0.697629772271099
		fp.875.n <- 1 - exp(-0.310122738776431 ) / n^ 1.06241615923172
	    }
	} else { ## --- p > 1 ---
	    if (intercept) {
		##			     "alfaq"		"betaq"	   "qwaarden"
		coefgqpkwad875 <- matrix(c(-0.458580153984614, 1.12236071104403, 3,
					   -0.267178168108996, 1.1022478781154,	 5), ncol = 2)
		coefeqpkwad500 <- matrix(c(-0.746945886714663, 0.56264937192689,  3,
					   -0.535478048924724, 0.543323462033445, 5), ncol = 2)
	    }
	    else {
		##			     "alfaq"		"betaq"	   "qwaarden"
		coefgqpkwad875 <- matrix(c(-0.251778730491252, 0.883966931611758, 3,
					   -0.146660023184295, 0.86292940340761,  5), ncol = 2)
		coefeqpkwad500 <- matrix(c(-0.487338281979106, 0.405511279418594, 3,
					   -0.340762058011,    0.37972360544988,  5), ncol = 2)
	    }

	    y.500 <- log(- coefeqpkwad500[1, ] / p^ coefeqpkwad500[2, ])
	    y.875 <- log(- coefgqpkwad875[1, ] / p^ coefgqpkwad875[2, ])

	    A.500 <- cbind(1, - log(coefeqpkwad500[3, ] * p^2))
	    coeffic.500 <- solve(A.500, y.500)
	    A.875 <- cbind(1, - log(coefgqpkwad875[3, ] * p^2))

	    coeffic.875 <- solve(A.875, y.875)
	    fp.500.n <- 1 - exp(coeffic.500[1]) / n^ coeffic.500[2]
	    fp.875.n <- 1 - exp(coeffic.875[1]) / n^ coeffic.875[2]
	}

	if(alpha <= 0.875)
	    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
	else ##	 0.875 < alpha <= 1
	    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)
    }## else (p >= 1)

    return(1/fp.alpha.n)
} ## LTScnp2

##' Compute Finite Sample Correction Factor for the  REWeighted LTSreg() scale
LTScnp2.rew <- function(p, intercept = intercept, n, alpha)
{
    stopifnot(0.5 <= alpha, alpha <= 1)
    if (intercept)
	p <- p - 1
    stopifnot(p == as.integer(p), p >= 0)

    if (p == 0) {
	fp.500.n <- 1 - exp( 1.11098143415027) / n^ 1.5182890270453
	fp.875.n <- 1 - exp(-0.66046776772861) / n^ 0.88939595831888

	if(alpha <= 0.875)
	    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
	else ##	 0.875 < alpha <= 1
	    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)
	## MM: sqrt() {below} is ''different logic'' than below.. (??)
	fp.alpha.n <- sqrt(fp.alpha.n)
    }
    else {
	if (p == 1) {
	    if (intercept) {
		fp.500.n <- 1 - exp(1.58609654199605 ) / n^ 1.46340162526468
		fp.875.n <- 1 - exp(0.391653958727332) / n^ 1.03167487483316
	    }
	    else {
		fp.500.n <- 1 - exp( 0.6329852387657)	/ n^ 1.40361879788014
		fp.875.n <- 1 - exp(-0.642240988645469) / n^ 0.926325452943084
	    }
	}
	else { ##  --- p > 1 ---
	    if (intercept) {
		##			     "alfaq"		"betaq"	   "qwaarden"
		coefqpkwad875 <- matrix(c(-0.474174840843602, 1.39681715704956, 3,
					  -0.276640353112907, 1.42543242287677, 5), ncol = 2)
		coefqpkwad500 <- matrix(c(-0.773365715932083, 2.02013996406346, 3,
					  -0.337571678986723, 2.02037467454833, 5), ncol = 2)
	    }
	    else {
		##			     "alfaq"		"betaq"	   "qwaarden"
		coefqpkwad875 <- matrix(c(-0.267522855927958, 1.17559984533974, 3,
					  -0.161200683014406, 1.21675019853961, 5), ncol = 2)
		coefqpkwad500 <- matrix(c(-0.417574780492848, 1.83958876341367, 3,
					  -0.175753709374146, 1.8313809497999, 5), ncol = 2)
	    }
	    y.500 <- log( - coefqpkwad500[1, ] / p^ coefqpkwad500[2, ])
	    y.875 <- log( - coefqpkwad875[1, ] / p^ coefqpkwad875[2, ])
	    A.500 <- cbind(1, - log(coefqpkwad500[3, ] * p^2))
	    coeffic.500 <- solve(A.500, y.500)
	    A.875 <- cbind(1, - log(coefqpkwad875[3, ] * p^2))
	    coeffic.875 <- solve(A.875, y.875)
	    fp.500.n <- 1 - exp(coeffic.500[1]) / n^ coeffic.500[2]
	    fp.875.n <- 1 - exp(coeffic.875[1]) / n^ coeffic.875[2]
	}

	if(alpha <= 0.875)
	    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
	else ##	 0.875 < alpha <= 1
	    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)

    }## else (p >= 1)

    return(1/fp.alpha.n)
} ## LTScnp2.rew

.fastlts <- function(x, y, h.alph, nsamp, intercept, adjust, trace = 0)
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    ## Parameters for partitioning --- *IDENTICAL* to those in ../src/rfltsreg.[fc]
    kmini <- 5
    nmini <- 300
    km10 <- 10*kmini
    nmaxi <- nmini*kmini

    ##	 vt::03.02.2006 - added options "best" and "exact" for nsamp
    if(!missing(nsamp)) {
	if(trace) cat("non-missing nsamp = ", nsamp, "\n")
	if(is.numeric(nsamp) && nsamp <= 0) {
	    warning("Invalid number of trials nsamp=",nsamp,"! Using default.\n")
	    nsamp <- -1
	} else if(nsamp == "exact" || nsamp == "best") {
	    myk <- p
	    if(n > 2*nmini-1) {
		warning("'nsamp' options 'best' and 'exact' not allowed for n greater than ",
                        2*nmini-1,". Will use default.\n")
		nsamp <- -1
	    }
            else { ## FIXME: Add a test case for this !
		nall <- choose(n, myk)
		if(nall > 5000 && nsamp == "best") {
		    nsamp <- 5000
		    warning("Maximum 5000 subsets allowed for option 'best'.\n",
                            "Computing 5000 subsets of size ",myk," out of ",n,"\n")
		} else {
		    nsamp <- 0		#all subsamples
		    if(nall > 5000)
			cat("Computing all ",nall," subsets of size ", myk,
                            " out of ",n,
                            "\n This may take a very long time!\n")
		}
	    }
        }
	if(nsamp == -1) { ## still not defined - set it to the default
	    nsamp <- rrcov.control()$nsamp
	}
    }
    nsamp <- as.integer(nsamp)

    ## y <- as.matrix(y)
    ## xy <- matrix(0, ncol = p + 1, nrow = n)
    xy <- cbind(x, y)
    storage.mode(xy) <- "double" # {keeping dim(.)}
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer" ; p1 <- p+1L # integer
    storage.mode(h.alph) <- "integer"

    ##	 Allocate temporary storage for the fortran implementation

    temp <- index1 <- index2 <- integer(n)

    weights <- aw2 <- aw <- residu <- yy <-
	nmahad <- ndist <- am <- am2 <- slutn <- double(n)

    .Fortran(rfltsreg, ## -> ../src/rfltsreg.f
	     xy = xy,
	     n,
	     p,
	     h.alph, # = nhalff
	     nsamp,  # = krep

	     inbest = integer(h.alph),
	     objfct = -1.,# double, if remains at -1 : have *nothing* found

	     intercept = as.integer(intercept),
	     intadjust = as.integer(adjust),
	     nvad = as.integer(p1),
	     datt = matrix(0., ncol = p1, nrow = n),
	     weights,
	     temp,
	     index1,
	     index2,
	     aw2,
	     aw,
	     residu,
	     yy,
	     nmahad,
	     ndist,
	     am, am2,
	     slutn,
             jmiss = integer(p1),	##	 integer jmiss(nvad)	  --> p+1
             xmed = double(p1),		##	 double	 xmed(nvad)	  --> p+1
             xmad = double(p1),		##	 double	 xmad(nvad)
             a	 = double(p1),		##	 double	    a(nvad)
             da	 = double(p1),		##	 double	   da(nvad)

             h = matrix(0., p, p1),	##	 double	 h(nvar,nvad)		p*(p+1)
             hvec = double(p*(p1)),	##	 double	 hvec(nvar*nvad)	p*(p+1)
             c = matrix(0., p, p1),	##	 double	 c(nvar,nvad)		p*(p+1)

             cstock = matrix(0., 10, p*p),##	 double	 cstock(10,nvar*nvar)	10*p*p
             mstock = matrix(0., 10, p), ##	 double	 mstock(10,nvar)	10*p
             c1stock =matrix(0., km10, p*p),##	 double	 c1stock(km10,nvar*nvar)  km10*p*p
             m1stock =matrix(0., km10, p),##	 double	 m1stock(km10,nvar)	km10*p

             dath  = matrix(0., nmaxi, p1),##	 double	 dath(nmaxi,nvad)	nmaxi*(p+1)
             sd    = double(p),		##	 double	 sd(nvar)		p
             means = double(p),		##	 double	 means(nvar)		p
             bmeans= double(p),		##	 double	 means(nvar)		p

         i.trace= as.integer(trace))[ c("inbest", "objfct") ]
}
