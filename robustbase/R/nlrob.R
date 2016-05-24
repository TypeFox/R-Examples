## Directly use nls()-internals, i.e., its 'm', to get a next 'start' (coef-like list):
## (In principle useful also outside robustbase)
.nls.get.start <- function(nls.m) {
    ## stopifnot(is.list(nls.m), is.function(gg <- nls.m$getPars),
    ##           is.environment(em <- environment(gg)))
    stopifnot(is.list(nls.m), is.environment(em <- environment(nls.m$getPars)))
    mget(names(em$ind), em$env)
}

nlrob <-
    function (formula, data, start, lower, upper,
              weights = NULL, na.action = na.fail,
	      method = c("M", "MM", "tau", "CM", "mtl"),
	      psi = .Mwgt.psi1("huber", cc=1.345), scale = NULL,
	      test.vec = c("resid", "coef", "w"),
	      maxit = 20, tol = 1e-06, acc,
	      algorithm = "default", doCov = FALSE,
	      control = if(method == "M") nls.control() else
			nlrob.control(method, optArgs = list(trace=trace), ...),
              trace = FALSE, ...)
{
    ## Purpose:
    ##	Robust fitting of nonlinear regression models. The fitting is
    ##	done by iterated reweighted least squares (IWLS) as in rlm() of
    ##	the package MASS. In addition, see also 'nls'.
    ##
    ## --> see the help file,  ?nlrob  (or ../man/nlrob.Rd in the source)
    ## -------------------------------------------------------------------------

    ##- some checks
    mf <- call <- match.call() # << and more as in nls()
    formula <- as.formula(formula)
    if (length(formula) != 3)
	stop("'formula' should be a formula of the type 'y ~ f(x, alpha)'")
    ## Had 'acc'; now use 'tol' which is more universal; 'acc' should work for a while
    if(!missing(acc) && is.numeric(acc)) {
        if(!missing(tol)) stop("specifying both 'acc' and 'tol' is invalid")
        tol <- acc
        message("The argument 'acc' has been renamed to 'tol'; do adapt your code.")
    }
    method <- match.arg(method)
    dataName <- substitute(data)
    dataCl <- attr(attr(call, "terms"), "dataClasses")
    hasWgts <- !missing(weights) # not eval()ing !
    if(method != "M") {
      if(hasWgts) ## FIXME .. should not be hard, e.g. for MM
          stop("specifying 'weights' is not yet supported for method ", method)
      if(!missing(psi))
	  warning(gettextf("For method = \"%s\", currently 'psi' must be specified via 'control'",
			   method), domain=NA)
      ## lifted from Martin's 'sfsmisc' package :
      missingCh <- function(x, envir = parent.frame()) {
          eval(substitute(missing(VAR), list(VAR=as.name(x))), envir = envir)
      }
      aNms <- c("start", "na.action", "test.vec", "maxit", "algorithm", "doCov")
      not.missA <- !vapply(aNms, missingCh, NA, envir=environment())
      if(any(not.missA)) {
	  warning(sprintf(ngettext(sum(not.missA),
				   "For method = \"%s\", argument %s is not made use of",
				   "For method = \"%s\", arguments %s are not made use of"),
			  method, pasteK(sQuote(aNms[not.missA]))),
		  domain=NA)
      }
      force(control)
      fixAns <- function(mod) {
          mod$call <- call # the nlrob() one, not nlrob.<foo>()
          ctrl <- mod$ctrl
          if(is.character(psi <- ctrl$psi) && is.numeric(cc <- ctrl$tuning.psi.M)) {# MM:
              psi <- .Mwgt.psi1(psi, cc=cc)
              res.sc <- with(mod, residuals/Scale)
              mod$psi <- psi
              mod$w <- # as we have no 'weights' yet
              mod$rweights <- psi(res.sc)
          } ## else mod$rweights <- mod$psi <- NULL
          mod$data <- dataName
          mod$dataClasses <- dataCl
          mod
      }
      switch(method,
	     "MM" = {
		 return(fixAns(nlrob.MM (formula, data, lower=lower, upper=upper,
					 tol=tol, ctrl= control)))
	     },
	     "tau" = {
		 return(fixAns(nlrob.tau(formula, data, lower=lower, upper=upper,
					 tol=tol, ctrl= control)))
	     },
	     "CM" = {
		 return(fixAns(nlrob.CM (formula, data, lower=lower, upper=upper,
					 tol=tol, ctrl= control)))
	     },
	     "mtl" = {
		 return(fixAns(nlrob.mtl(formula, data, lower=lower, upper=upper,
					 tol=tol, ctrl= control)))
	     })
    } ## {non-"M" methods}
    ## else: method == "M", original method, the only one based on 'nls' :

    varNames <- all.vars(formula)
    env <- environment(formula)
    if (is.null(env)) env <- parent.frame()
    ## FIXME:  nls() allows  a missing 'start'; we don't :
    if(length(pnames <- names(start)) != length(start))
        stop("'start' must be fully named (list or numeric vector)")
    if (!((is.list(start) && all(sapply(start, is.numeric))) ||
	  (is.vector(start) && is.numeric(start))))
	stop("'start' must be a named list or numeric vector")
    if(any(is.na(match(pnames, varNames))))
	stop("parameter names must appear in 'formula'")
    ## If it is a parameter it is not a variable
    varNames <- varNames[is.na(match(varNames, pnames))]
    test.vec <- match.arg(test.vec)
    if(missing(lower)) lower <- -Inf
    if(missing(upper)) upper <- +Inf
    updateScale <- is.null(scale)
    if(!updateScale) { ## keep initial scale fixed through iterations (e.g. for "MM")
	if(is.1num(scale) && scale > 0)
            Scale <- scale
        else
            stop("'scale' must be NULL or a positive number")
    }
    nm <- "._nlrob.w"
    if (nm %in% c(varNames, pnames, names(data)))
	stop(gettextf("Do not use '%s' as a variable name or as a parameter name",
		      nm), domain=NA)

    obsNames <- rownames(data <- as.data.frame(data))
    data <- as.list(data)# to be used as such

    ## From nls: using model.weights() e.g. when formula 'weights = sqrt(<var>)'
    mf$formula <-  # replace by one-sided linear model formula
        as.formula(paste("~", paste(varNames, collapse = "+")),
                   env = environment(formula))
    mf[c("start", "lower", "upper", "method", "psi", "scale", "test.vec",
         "maxit", "tol", "acc", "algorithm", "doCov", "control", "trace")] <- NULL
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(mf)
    nobs <- nrow(mf)
    mf <- as.list(mf)
    if (hasWgts)
	hasWgts <- !is.null(weights <- model.weights(mf))
    if (hasWgts && any(weights < 0 | is.na(weights)))
	stop("'weights' must be nonnegative and not contain NAs")
    ## initialize testvec etc
    fit <- eval(formula[[3]], c(data, start), env)
    y <- eval(formula[[2]], data, env)
    coef <- unlist(start)
    resid <- y - fit
    ## if (any(is.na(data)) & options("na.action")$na.action == "na.omit")
    ##	 stop("if NAs are present, use 'na.exclude' to preserve the residuals length")

    irls.delta <- function(old, new) sqrt(sum((old - new)^2, na.rm = TRUE)/
					  max(1e-20, sum(old^2, na.rm = TRUE)))
    ## Robust loop -- IWLS / IRLS iterations
    converged <- FALSE
    status <- "converged"
    method.exit <- FALSE
    for (iiter in seq_len(maxit)) {
	if (trace)
	    cat("robust iteration", iiter, "\n")
	previous <- get(test.vec)
	if(updateScale)
            Scale <- median(abs(resid), na.rm = TRUE)/0.6745
	if (Scale == 0) {
	    convi <- 0
	    method.exit <- TRUE
	    warning(status <- "could not compute scale of residuals")
	    ## FIXME : rather use a "better" Scale in this case, e.g.,
	    ## -----   Scale <- min(abs(resid)[resid != 0])
	}
	else {
	    w <- psi(resid/Scale)
	    if (hasWgts)
		w <- w * weights
	    data$._nlrob.w <- w ## use a variable name the user "will not" use
	    ._nlrob.w <- NULL # workaround for codetools "bug"
            ## Case distinction against "wrong warning" as long as
            ## we don't require R > 3.0.2:
            out <-
                if(identical(lower, -Inf) && identical(upper, Inf))
                    nls(formula, data = data, start = start,
                        algorithm = algorithm, trace = trace,
                        weights = ._nlrob.w,
                        na.action = na.action, control = control)
                else
                    nls(formula, data = data, start = start,
                        algorithm = algorithm, trace = trace,
                        lower=lower, upper=upper,
                        weights = ._nlrob.w,
                        na.action = na.action, control = control)

            coef <- unlist(start <- .nls.get.start(out$m))
	    ## same sequence as in start! Ok for test.vec:
            resid <- residuals(out)
	    convi <- irls.delta(previous, get(test.vec))
	}
	converged <- convi <= tol
	if (converged)
	    break
	else if (trace)
	    cat(sprintf(" --> irls.delta(previous, %s) = %g -- *not* converged\n",
                        test.vec, convi))
    }## for( iiter ...)

    if(!converged || method.exit) {
	warning(st <- paste("failed to converge in", maxit, "steps"))
        status <- if(method.exit) {
            converged <- FALSE; paste(status, st, sep="; ") } else st
    }

    if(hasWgts) { ## or just   out$weights  ??
	tmp <- weights != 0
	w[tmp] <- w[tmp]/weights[tmp]
    }

    ## --- Estimated asymptotic covariance of the robust estimator
    rw <- psi(res.sc <- resid/Scale)
    asCov <- if(!converged || !doCov) NA else {
        ## a version of  .vcov.m(.) below
	AtWAinv <- chol2inv(out$m$Rmat())
	dimnames(AtWAinv) <- list(names(coef), names(coef))
	tau <- mean(rw^2) / mean(psi(res.sc, d=TRUE))^2
	AtWAinv * Scale^2 * tau
    }

    ## returned object:	 ==  out$m$fitted()  [FIXME?]
    fit <- setNames(eval(formula[[3]], c(data, start)), obsNames)
    structure(class = c("nlrob", "nls"),
	      list(m = out$m, call = call, formula = formula,
		   new.formula = formula, nobs = nobs,
		   coefficients = coef,
		   working.residuals = as.vector(resid),
		   fitted.values = fit, residuals = y - fit,
		   Scale=Scale, w=w, rweights = rw,
		   cov = asCov, test.vec=test.vec, status=status, iter=iiter,
		   psi=psi, data = dataName, dataClasses = dataCl,
		   control = control))
}

##' @title The nlrob() method used
##' @param obj an \code{"nlrob"} object
##' @return characer string
.method.nlrob <- function(obj) if(inherits(obj, "nls")) "M" else obj$ctrl$method

.vcov.m <- function(object, Scale, resid.sc) {
    if(.method.nlrob(object) == "M") {
	AtWAinv <- chol2inv(object$m$Rmat())
	stopifnot(length(Scale) == 1, Scale >= 0,
		  is.numeric(resid.sc), length(resid.sc) == nobs(object),
		  is.character(nms.coef <- names(coef(object))),
		  length(nms.coef) == nrow(AtWAinv),
		  is.function(psi <- object$psi))

	dimnames(AtWAinv) <- list(nms.coef, nms.coef)
	tau <- mean(psi(resid.sc)^2) / mean(psi(resid.sc, d=TRUE))^2
	AtWAinv * Scale^2 * tau
    }
    else if(is.function(psi <- object$psi)) {
	form <- object$formula
	## call method="M", with fixed Scale
	mM <- nlrob(form, data = eval(object$data, environment(form)),
		    method = "M", start = coef(object),
		    psi = psi, scale = Scale, doCov=TRUE)
	mM$cov
	## stop(".vcov.m() not yet implemented for nlrob.MM objects")
	## using 'chol(<Hessian>): --- is wrong, unfortunately
	## AtWAinv <- chol2inv(chol(object$hessian))
    } else {
	NA ## instead of error
    }
}


## The 'nls' method is *not* correct
formula.nlrob <- function(x, ...) x$formula

sigma.nlrob <- function(object, ...)
    if(!is.null(s <- object$Scale)) s else object$coefficients[["sigma"]]

estimethod <- function(object, ...) UseMethod("estimethod")
estimethod.nlrob <- function(object, ...)
    if(is.list(object$m) && inherits(object, "nls")) "M" else object$ctrl$method

fitted.nlrob <- function (object, ...)
{
    val <- as.vector(object$fitted.values)
    if (!is.null(object$na.action))
	val <- napredict(object$na.action, val)
    ##MM: attr(val, "label") <- "Fitted values"
    val
}


## formula() works "by default"

predict.nlrob <- function (object, newdata, ...)
{
    if (missing(newdata))
	return(as.vector(fitted(object)))
    if (!is.null(cl <- object$dataClasses))
	.checkMFClasses(cl, newdata)
    if(estimethod(object) == "M") # also for start = list(..)
	object$m$predict(newdata)
    else
	eval(formula(object)[[3]], c(as.list(newdata), coef(object)))
}


print.nlrob <- function (x, ...)
{
    cat("Robustly fitted nonlinear regression model",
	if((meth <- .method.nlrob(x)) != "M") paste0(" (method ", meth, ")"),
	"\n", sep="")
    cat("  model: ", deparse(formula(x)), "\n")
    cat("   data: ", deparse(x$data), "\n")
    print(coef(x), ...)
    cat(" status: ", x$status, "\n")
    invisible(x)
}


residuals.nlrob <- function (object, type = c("response", "working", "pearson"), ...)
{
    type <- match.arg(type)
    R <- switch(type,
                "pearson"=
            {
                stop("type 'pearson' is not yet implemented")
                ## as.vector(object$working.residuals)
            },
                "working"=
            {   ## FIXME(?): from nls, these used to *contain* weights, but no longer
                object$working.residuals
            },
                "response"=
            {
                object$residuals
            },
                stop("invalid 'type'"))# ==> programming error, as we use match.arg()
    if (!is.null(object$na.action))
        R <- naresid(object$na.action, R)
    ## FIXME: add 'names'!
    ##MM no labels; residuals.glm() does neither: attr(val, "label") <- "Residuals"
    R
}


vcov.nlrob <- function (object, ...) {
    if(is.numeric(cv <- object$cov)) cv
    else {
        sc <- object$Scale
        .vcov.m(object, Scale = sc, resid.sc = as.vector(object$residuals) / sc)
    }
}

summary.nlrob <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    w <- object$w ## weights * rweights, scaled such that sum(w)=1
    n <- sum(w > 0)
    param <- coef(object)
    p <- length(param)
    rdf <- n - p
    no <- names(object)
    no <- no[match(c("formula", "residuals", "Scale", "w", "rweights", "cov",
                     "call", "status", "counts", "iter", "control", "ctrl"), no, 0L)]
    ans <- object[no]
    conv <- ans$status == "converged"
    sc   <- ans$Scale
    if(conv && !is.matrix(ans$cov))
	ans$cov <- .vcov.m(object, Scale = sc,
			   resid.sc = as.vector(object$residuals) / sc)
    if((ok.cov <- is.matrix(ans$cov)))
        if(!all(dim(ans$cov) == p)) stop("'cov' must be a p x p matrix")
    ans$df <- c(p, rdf)
    cf <-
	if(ok.cov) {
	    se <- sqrt(diag(ans$cov))
	    tval <- param/se
	    cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
	} else cbind(param, NA, NA, NA)
    dimnames(cf) <- list(names(param),
			 c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$coefficients <- cf
    if(correlation && ok.cov && rdf > 0) {
	ans$correlation <- ans$cov / outer(se, se)
	ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.nlrob"
    ans
}

print.summary.nlrob <-
    function (x, digits = max(3, getOption("digits") - 3),
            symbolic.cor = x$symbolic.cor,
            signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
	"\n\n", sep = "")
    ## cat("\nFormula: ")
    ## cat(paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n", sep = "")
    if(is.null(ctrl <- x$ctrl))
        meth <- "M"
    else {
	meth <- ctrl$method
	cat("Method \"", meth,
	    if(!is.null(cc <- ctrl$init)) paste0("\", init = \"", cc),
	    if(!is.null(ps <- ctrl$psi )) paste0("\", psi = \"", ps),
	    "\"\n", sep="")
    }
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
    cat("\nParameters:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
		 ...)
    if(x$status == "converged") {
	cat("\nRobust residual standard error:",
	    format(signif(x$Scale, digits)), "\n")
	correl <- x$correlation
	if (!is.null(correl)) {
	    p <- NCOL(correl)
	    if (p > 1) {
		cat("\nCorrelation of Parameter Estimates:\n")
		if(is.logical(symbolic.cor) && symbolic.cor) {
		    print(symnum(correl, abbr.colnames = NULL))
		} else {
		    correl <- format(round(correl, 2), nsmall = 2, digits = digits)
		    correl[!lower.tri(correl)] <- ""
		    print(correl[-1, -p, drop=FALSE], quote = FALSE)
		}
	    }
	}
	if(is.null(ctrl))
	    cat("Convergence in", x$iter, "IRWLS iterations\n\n")
	else {
	    if(length(it <- ctrl$iter) == 1)
		cat("Convergence in", it, "iterations\n\n")
	    else if(length(cnts <- x$counts) > 0)
		cat("Convergence after", cnts[["function"]],
		    "function and", cnts[["gradient"]],"gradient evaluations\n\n")
	    else ## length(it) >= 2 :
		cat("Convergence\n\n")
	}
	summarizeRobWeights(x$rweights, digits = digits, ...)
    }
    else if(meth == "M")
	cat("** IRWLS iterations did *not* converge!\n\n")
    else
	cat("** Iterations did *not* converge!\n\n")
    invisible(x)
}

## Confint(): ideally built on profile, the same as stats:::confint.nls()
## --------   which eventually calls stats:::profile.nls

## Also, do emulate (to some extent)

## str(lme4:::confint.merMod)
## function (object, parm, level = 0.95, method = c("profile", "Wald", "boot"),
##     zeta, nsim = 500, boot.type = c("perc", "basic", "norm"), quiet = FALSE,
##     oldNames = TRUE, ...)

confint.nlrob <-
    function(object, parm, level = 0.95,
	     method = c("profile", "Wald", "boot"),
	     zeta, nsim = 500, boot.type = c("perc", "basic", "norm"),
	     quiet = FALSE, oldNames = TRUE, ...)
{
    method <- match.arg(method)
    boot.type <- match.arg(boot.type)
    if (!missing(parm) && !is.numeric(parm) &&
	method %in% c("profile", "boot"))
	stop("for method='", method, "', 'parm' must be specified as an integer")
    switch(method, profile = {
	stop("profile() method not yet implemented for \"nlrob\" objects.
 Use  method = \"Wald\".")
	## hence unused for now :
	if (!quiet) message("Computing profile confidence intervals ...")
	utils::flush.console()
	pp <- if (missing(parm)) {
	    profile(object, signames = oldNames, ...)
	} else {
	    profile(object, which = parm, signames = oldNames,
		...)
	}
	confint(pp, level = level, zeta = zeta)
    }, Wald = {
	cf <- coef(object)
	pnames <- names(cf)
	if (missing(parm)) parm <- pnames
	else if (is.numeric(parm)) parm <- pnames[parm]
	a <- (1 - level)/2
	a <- c(a, 1 - a)
	## for now, a short version of R's formatting in quantile.default():
	format_perc <- function(x, digits = max(2L, getOption("digits")))
	    paste0(formatC(x, format = "fg", width = 1, digits = digits))
	pct <- format_perc(a, 3)
	fac <- qnorm(a)
	ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
	sdiag <- function(x) if (length(x) == 1) c(x) else diag(x)
	ses <- sqrt(sdiag(vcov(object)[parm, parm]))
	ci[] <- cf[parm] + ses %o% fac
	ci
    }, boot = {
	stop("\"boot\" method not yet implemented for \"nlrob\" objects.
 Use confint(*, method = \"Wald\").")
    })
}
