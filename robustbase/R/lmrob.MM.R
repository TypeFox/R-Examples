## The "regularized" psi-function names:
## .R: the redescending ones:
.Mpsi.R.names <- c('bisquare', 'lqq', 'welsh', 'optimal', 'hampel', 'ggw')
## .M: the monotone ones:
.Mpsi.M.names <- c('huber')
## Note: there could be more: non-redescending, non-monotone {such as Cauchy score}
.Mpsi.names <- c(R= .Mpsi.R.names, M= .Mpsi.M.names)


##' This allows synonyms as "Tukey" *and* partial matches such as "opt" :
.regularize.Mpsi <- function(psi, redescending = TRUE) {
    stopifnot(is.character(psi), length(psi) == 1)
    psi <- tolower(psi)
    psi <- switch(psi,
		  'tukey'= , 'biweight'= "bisquare",
		  ## otherwise keep
		  psi)
    nms <- if(redescending) .Mpsi.R.names else .Mpsi.names
    if (is.na(i <- pmatch(psi, nms)))
	stop(gettextf("'psi' should be one of %s", pasteK(dQuote(nms))),
	     domain = NA)
    nms[i]
}

.Mpsi.tuning.defaults <- list(
    'huber' = 1.345
    , 'bisquare' = 4.685061
    , 'welsh' = 2.11
    , 'ggw' = c(-0.5, 1.5, .95, NA) ## min slope, b, eff, bp
    , 'lqq' = c(-0.5, 1.5, .95, NA) ## min slope, b, eff, bp
    , 'optimal' = 1.060158
    , 'hampel' = c(1.5, 3.5, 8) * 0.9016085 ## a, b, r
    )
.Mpsi.tuning.default <- function(psi) {
    if(is.null(p <- .Mpsi.tuning.defaults[[psi]]))
        stop(gettextf("invalid 'psi'=%s; possibly use .regularize.Mpsi(%s)",
                      psi, "psi, redescending=FALSE"), domain=NA)
    p
}

.Mchi.tuning.defaults <- list(
    ## Here, psi must be redescending! -> 'huber' not possible
    'bisquare' = 1.54764
    , 'welsh' = 0.5773502
    , 'ggw' = c(-0.5, 1.5, NA, .5) ## min slope, b, eff, bp
    , 'lqq' = c(-0.5, 1.5, NA, .5) ## min slope, b, eff, bp
    , 'optimal' = 0.4047
    , 'hampel' = c(1.5, 3.5, 8) * 0.2119163 ## a, b, r
    )
.Mchi.tuning.default <- function(psi) {
    if(is.null(p <- .Mchi.tuning.defaults[[psi]]))
	stop(gettextf("invalid 'psi'=%s; possibly use .regularize.Mpsi(%s)",
		      psi, "psi"), domain=NA)
    p
}


lmrob.control <-
    function(setting, seed = NULL, nResample = 500,
	     tuning.chi = NULL, bb = 0.5,
	     tuning.psi = NULL, max.it = 50,
	     groups = 5, n.group = 400, k.fast.s = 1, best.r.s = 2,
	     k.max = 200, maxit.scale = 200, k.m_s = 20,
	     ##           ^^^^^^^^^^^ had MAX_ITER_FIND_SCALE 200 in ../src/lmrob.c
	     refine.tol = 1e-7, rel.tol = 1e-7,
	     solve.tol = 1e-7,
	     ## had  ^^^^^^^^  TOL_INVERSE 1e-7 in ../src/lmrob.c
	     trace.lev = 0, mts = 1000,
	     subsampling = c("nonsingular", "simple"),
	     compute.rd = FALSE,
	     method = 'MM',
	     psi = 'bisquare',
	     numpoints = 10, cov = NULL,
	     split.type = c("f", "fi", "fii"),
	     fast.s.large.n = 2000,
             eps.outlier = function(nobs) 0.1 / nobs,
             eps.x = function(maxx) .Machine$double.eps^(.75)*maxx,
             compute.outlier.stats = method,
             warn.limit.reject = 0.5,
             warn.limit.meanrw = 0.5,
             ...)
{
    p.ok <- missing(psi) # if(p.ok) psi does not need regularization
    if (!missing(setting)) {
        if (setting %in% c('KS2011', 'KS2014')) {
            if (missing(method)) method <- 'SMDM'
	    psi <- if(p.ok) 'lqq' else .regularize.Mpsi(psi) ; p.ok <- TRUE
            if (missing(max.it)) max.it <- 500
            if (missing(k.max)) k.max <- 2000
            if (missing(cov) || is.null(cov)) cov <- '.vcov.w'
            if (setting == 'KS2014') {
                if (missing(best.r.s)) best.r.s <- 20
                if (missing(k.fast.s)) k.fast.s <- 2
                if (missing(nResample)) nResample <- 1000
            }
        } else {
            warning("Unknown setting '", setting, "'. Using defaults.")
        }
    } else {
	if(p.ok && grepl('D', method)) psi <- 'lqq'
	if (missing(cov) || is.null(cov))
	    cov <- if(method %in% c('SM', 'MM')) ".vcov.avar1" else ".vcov.w"
    }
    if(!p.ok) psi <- .regularize.Mpsi(psi)
    subsampling <- match.arg(subsampling)

    ## in ggw, lqq:  if tuning.{psi|chi}  are non-standard, calculate coefficients:
    compute.const <- (psi %in% c('ggw', 'lqq'))

    if(is.null(tuning.chi))
	tuning.chi <- .Mchi.tuning.default(psi)
    else if(compute.const)
	tuning.chi <- .psi.const(tuning.chi, psi)

    if(is.null(tuning.psi))
	tuning.psi <- .Mpsi.tuning.default(psi)
    else if(compute.const)
	tuning.psi <- .psi.const(tuning.psi, psi)

    c(list(setting = if (missing(setting)) NULL else setting,
           seed = as.integer(seed), nResample=nResample, psi=psi,
           tuning.chi=tuning.chi, bb=bb, tuning.psi=tuning.psi,
           max.it=max.it, groups=groups, n.group=n.group,
           best.r.s=best.r.s, k.fast.s=k.fast.s,
           k.max=k.max, maxit.scale=maxit.scale, k.m_s=k.m_s, refine.tol=refine.tol,
           rel.tol=rel.tol, solve.tol=solve.tol, trace.lev=trace.lev, mts=mts,
           subsampling=subsampling,
           compute.rd=compute.rd, method=method, numpoints=numpoints,
           cov=cov, split.type = match.arg(split.type),
           fast.s.large.n=fast.s.large.n,
           eps.outlier = eps.outlier, eps.x = eps.x,
           compute.outlier.stats = sub("^MM$", "SM", compute.outlier.stats),
           warn.limit.reject = warn.limit.reject,
           warn.limit.meanrw = warn.limit.meanrw),
      list(...))
}

##' Modify a \code{\link{lmrob.control}} list to contain only parameters that
##' were actually used.  Currently used for \code{\link{print}()}ing of lmrob
##' objects.
##'
##' @title Minimize lmrob control to non-redundant parts
##' @param control a list, typically the 'control' component of a
##' \code{\link{lmrob}()} call, or the result of  \code{\link{lmrob.control}()}.
##' @return list: the (typically) modified \code{control}
##' @author Martin Maechler {from Manuel's original code}
lmrob.control.neededOnly <- function(control) {
    if(is.null(control)) return(control)
    switch(sub("^(S|M-S).*", "\\1", control$method),
	   S = {                       # remove all M-S specific control pars
	       control$k.m_s <- NULL
	       control$split.type <- NULL
					# if large_n is not used, remove corresp control pars
	       if (length(residuals) <= control$fast.s.large.n) {
		   control$groups <- NULL
		   control$n.group <- NULL
	       }
	   },
	   `M-S` = {                # remove all fast S specific control pars
	       control$refine.tol <- NULL
	       control$groups <- NULL
	       control$n.group <- NULL
	       control$best.r.s <- NULL
	       control$k.fast.s <- NULL
	   }, { # else: do not keep parameters used by initial ests. only
	       control$tuning.chi <- NULL
	       control$bb <- NULL
	       control$refine.tol <- NULL
	       control$nResample <- NULL
	       control$groups <- NULL
	       control$n.group <- NULL
	       control$best.r.s <- NULL
	       control$k.fast.s <- NULL
	       control$k.max <- NULL
	       control$k.m_s <- NULL
	       control$split.type <- NULL
	       control$mts <- NULL
	       control$subsampling <- NULL
	   } )
    if (!grepl("D", control$method)) control$numpoints <- NULL
    if (control$method == 'SM') control$method <- 'MM'
    control
}

lmrob.fit.MM <- function(x, y, control) ## deprecated
{
    .Deprecated("lmrob.fit(*, control) with control$method = 'SM'")
    control$method <- 'SM'
    lmrob.fit(x, y, control)
}## lmrob.MM.fit()

lmrob.fit <- function(x, y, control, init=NULL, mf=NULL) {
    if(!is.matrix(x)) x <- as.matrix(x)
    ## old notation: MM -> SM
    if (control$method == "MM") control$method <- "SM"
    ## Assumption:  if(is.null(init))  method = "S..."   else  method = "..."
    ## ---------    where "..." consists of letters {"M", "D"}
    est <- if (is.null(init)) {
        ## --- initial S estimator
        if ((M1 <- substr(control$method,1,1)) != 'S') {
	    warning(gettextf("Initial estimator '%s' not supported; using S-estimator instead",
			     M1), domain = NA)
            substr(control$method,1,1) <- 'S'
        }
        init <- lmrob.S(x, y, control = control, mf = mf)
        'S'
    } else {
	stopifnot(is.list(init))
        if (is.null(init$converged)) init$converged <- TRUE
	if (is.null(init$control)) {
	    init$control <- control
	    M <- init$control$method <- 'l'
	} else if(!length(M <- init$control$method) || !nzchar(M))
	    M <- "l"
	M
    }
    stopifnot(is.numeric(init$coef), length(init$coef) == ncol(x),
              is.numeric(init$scale), init$scale >= 0)
    if (est != 'S' && control$cov == '.vcov.avar1') {
        warning(
	".vcov.avar1 can only be used when initial estimator is S; using .vcov.w instead")
        control$cov <- ".vcov.w"
    }
    if (init$converged) {
        ## --- loop through the other estimators; build up 'est' string
        method <- sub(paste0("^", est), '', control$method)
        for (step in strsplit(method,'')[[1]]) {
            ## now we have either M or D steps
            est <- paste0(est, step)
            init <- switch(step,
                           ## D(AS)-Step
                           D = lmrob..D..fit(init, x, mf = mf),
                           ## M-Step
                           M = lmrob..M..fit(x = x, y = y, obj = init, mf = mf),
                           stop('only M and D are steps supported after "init" computation'))
            ## break if an estimator did not converge
            if (!init$converged) {
		warning(gettextf(
		    "%s-step did NOT converge. Returning unconverged %s-estimate",
		    step, est),
			domain = NA)
                break
            }
        }
    }
    ## << FIXME? qr(.)  should be available from earlier
    if (is.null(init$qr)) init$qr <- qr(x * sqrt(init$rweights))
    if (is.null(init$rank)) init$rank <- init$qr$rank
    control$method <- est ## ~= original 'method', but only with the steps executed.
    init$control <- control

    ## --- covariance estimate
    init$cov <-
	if (init$scale == 0) { ## exact fit
	    matrix(0, ncol(x), ncol(x),
		   dimnames=list(colnames(x), colnames(x)))
	} else if (!init$converged || is.null(x)) {
	    NA
	} else {
	    if (is.null(control$cov) || control$cov == "none")
		NA
	    else {
		lf.cov <- if (!is.function(control$cov))
		    get(control$cov, mode='function') else control$cov
		lf.cov(init, x)
	    }
	}
    df <- NROW(y) - init$rank ## sum(init$r?weights)-init$rank
    init$degree.freedom <- init$df.residual <- df
    init
}

globalVariables("r", add=TRUE) ## below and in other lmrob.E() expressions

.vcov.w <- function(obj, x=obj$x, scale=obj$scale, cov.hubercorr=ctrl$cov.hubercorr,
                    cov.dfcorr=ctrl$cov.dfcorr, cov.resid=ctrl$cov.resid,
                    cov.corrfact=ctrl$cov.corrfact,
                    cov.xwx=ctrl$cov.xwx)
{
    ## set defaults
    ctrl <- obj$control
    if (is.null(cov.hubercorr)) cov.hubercorr <- !grepl('D', ctrl$method)
    else if (!is.logical(cov.hubercorr))
        stop(':.vcov.w: cov.hubercorr must be logical (or NULL)')
    val.corrf <- c('tau', 'empirical', 'asympt', 'hybrid', 'tauold')
    if (is.null(cov.corrfact)) {
        cov.corrfact <- if (cov.hubercorr) 'empirical' else 'tau'
    } else if(length(cov.corrfact) != 1 || is.na(match(cov.corrfact, val.corrf)))
	stop(":.vcov.w: cov.corrfact must be one of ",
             pasteK(dQuote(val.corrf)))
    if (is.null(cov.dfcorr)) {
        cov.dfcorr <- if (cov.hubercorr | cov.corrfact %in% c('tau', 'hybrid')) 1 else -1
    } else if (!is.numeric(cov.dfcorr) || is.na(match(cov.dfcorr, -1:3)))
        stop(':.vcov.w: cov.dfcorr has to be one of -1:3')

    val.cov.res <- c('final', 'initial', 'trick')
    if (is.null(cov.resid)) cov.resid <- 'final'
    else if (length(cov.resid) != 1 || is.na(match(cov.resid, val.cov.res)))
	stop(":.vcov.w: cov.resid must be one of ",
	     pasteK(dQuote(val.cov.res)))
    if (is.null(cov.xwx)) cov.xwx <- TRUE
    else if (!is.logical(cov.xwx))
	stop(':.vcov.w: cov.xwx must be logical (or NULL)')
    if (is.null(x))  x <- model.matrix(obj)
    ## set psi and c.psi
    if (cov.resid == 'initial') {
        psi <- ctrl$psi
        c.psi <- ctrl$tuning.chi
        if (is.null(psi)) stop('parameter psi is not defined')
        if (!is.numeric(c.psi)) stop('parameter tuning.psi is not numeric')
    } else {
        psi <- ctrl$psi
        c.psi <- if (ctrl$method %in% c('S', 'SD'))
            ctrl$tuning.chi else ctrl$tuning.psi
        if (is.null(psi)) stop('parameter psi is not defined')
        if (!is.numeric(c.psi)) stop('parameter tuning.psi is not numeric')
    }
    if (cov.resid == 'final' && (class(obj)[1] == 'lmrob.S'))
        warning(":.vcov.w: ignoring cov.resid == final since est != final")
    if (is.null(scale)) {
        warning(":.vcov.w: scale missing, using D scale")
        scale <- lmrob..D..fit(obj)$scale
    }
    n <- NROW(x)
    ## --- calculations: matrix part
    ## weighted xtx.inv matrix
    w <- if (cov.xwx) obj$rweights else rep(1,n)
    ## use qr-decomposition from lm.wfit (this already includes the robustness weights)
    ## update qr decomposition if it is missing or we don't want the robustness weights
    if (!is.qr(obj$qr) || !cov.xwx) obj$qr <- qr(x * sqrt(w))
    p <- if (is.null(obj$rank)) obj$qr$rank else obj$rank
    cinv <- if(is.qr(obj$qr)) tryCatch(tcrossprod(solve(qr.R(obj$qr))),
				       error = function(e)e)
    if(inherits(cinv, 'error')) cinv <- matrix(NA,p,p)
    ## --- calculation: correction factor
    if (cov.corrfact == 'asympt') { ## asympt correction factor
        if (cov.hubercorr) warning('option hcorr is ignored for cov.corrfact = asympt')
        ## precalculated default values if applicable
        corrfact <-
            if (psi == 'ggw') {
                if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.95, NA)))) 1.052619
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.95, NA)))) 1.0525888644
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.85, NA)))) 1.176479
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.85, NA)))) 1.176464
                else lmrob.E(psi(r)^2, ctrl) / lmrob.E(r*psi(r), ctrl)^2
                ## MK: using r*psi(r) instead of psi'(r) is much more accurate
                ##     when using Gauss-Hermite quadrature
            } else if (isTRUE(all.equal(c.psi, .Mpsi.tuning.default(psi)))) {
                switch(psi,
                       bisquare = 1.0526317574,
                       welsh    = 1.0526704649,
                       optimal  = 1.0526419204,
                       hampel   = 1.0526016980,
                       lqq      = 1.0526365291,
                       stop(':.vcov.w: unsupported psi function'))
            } else lmrob.E(psi(r)^2, ctrl) / lmrob.E(r*psi(r), ctrl)^2 ## see above
        varcorr <- 1
    } else { ## empirical, approx or hybrid correction factor
	rstand <- if (cov.resid == 'initial') {
	    ## if the last estimator was a D or T estimator
	    ## then use obj$init$init otherwise use obj$init
	    ## that way for SMD we use the S residuals (and S scale)
	    ## and for SMDM we use the M residuals (and D scale)
	    lobj <-
		if (grepl('[DT]$',ctrl$method)) obj$init$init else obj$init
	    resid(lobj) / lobj$scale
	} else if (cov.resid == 'trick') {
	    ## residuals are in fact from earlier estimator, use its scale to standardize them
	    obj$init$resid / obj$init$scale
	} else obj$resid / scale

        tau <- if (cov.corrfact %in% c('tau', 'hybrid', 'tauold')) { ## added hybrid here
            if (!is.null(obj$tau)) obj$tau
            else if (!is.null(obj$init$tau)) obj$init$tau
            else stop("(tau / hybrid / tauold): tau not found in 'obj'") } else rep(1,n)
	rstand <- rstand / tau
        r.psi   <- Mpsi(rstand, c.psi, psi)
        r.psipr <- Mpsi(rstand, c.psi, psi, deriv = 1)
        if (any(is.na(r.psipr))) warning(":.vcov.w: Caution. Some psiprime are NA")
        mpp2 <- (mpp <- mean(r.psipr, na.rm=TRUE))^2
        ## Huber's correction
        hcorr <- if (cov.hubercorr) {
            vpp <- sum((r.psipr - mpp)^2) / n # var(r.psipr, na.rm=TRUE)
            (1+p/n*vpp/mpp2)^2
        } else 1
        ## sample size correction for var(r.psi^2)
        ## use tau if 'tau' correction factor, but only if it is available
        varcorr <- if (cov.corrfact == 'tau' && any(tau != 1))
            1 / mean(tau^2) else n / (n - p) ## changed from 1 / mean(tau)
        ## if hybrid: replace B (= mpp2) by asymptotic value
        if (cov.corrfact == 'hybrid') {
            mpp2 <- if (psi == 'ggw') {
                if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.95, NA)))) 0.7598857
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.95, NA)))) 0.6817983
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.85, NA)))) 0.4811596
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.85, NA)))) 0.411581
                else lmrob.E(r*psi(r), ctrl)^2 ## more accurate than psi'(r)
            } else if (isTRUE(all.equal(c.psi, lmrob.control(psi = psi)$tuning.psi)))
                switch(psi,
                       bisquare = 0.5742327, welsh = 0.5445068, optimal = 0.8598825,
                       hampel = 0.6775217, lqq = 0.6883393,
                       stop(':.vcov.w: unsupported psi function'))
            else lmrob.E(r*psi(r), ctrl)^2 ## more accurate than psi'(r)
        }
        corrfact <- mean({ if (cov.corrfact == 'tauold') 1 else tau^2 } * r.psi^2)/mpp2 * hcorr
    }
    ## simple sample size correction
    sscorr <- if (cov.dfcorr > 0) {
        if (cov.dfcorr == 2) varcorr ## cov.dfcorr == 2
        else if (cov.dfcorr == 3) mean(w)^2 / (1 - p / sum(w)) ## cov.dfcorr == 3
        else mean(w) * varcorr ## cov.dfcorr == 1
    } else if (cov.dfcorr < 0) mean(w) ## cov.dfcorr == -1
    else 1 ## cov.dfcorr == 0

    ## scale^2 * a/b2 * Huber's correction * Cinv
    cv <- scale^2 * sscorr * corrfact * cinv
    attr(cv,"weights") <- w
    attr(cv,"scale") <- scale
    attr(cv,"scorr") <- sscorr
    attr(cv,"corrfact") <- corrfact
    cv
}

.vcov.avar1 <- function(obj, x=obj$x, posdef.meth = c("posdefify","orig"))
{ ## was .vcov.MM
    stopifnot(is.list(ctrl <- obj$control))
    ## works only for MM & SM estimates:
    if (!is.null(ctrl$method) && !ctrl$method %in% c('SM', 'MM'))
        stop('.vcov.avar1() supports only SM or MM estimates')
    ## set psi and chi constants
    psi <- chi <- ctrl$psi
    if (is.null(psi)) stop('parameter psi is not defined')
    stopifnot(is.numeric(c.chi <- ctrl$tuning.chi),
	      is.numeric(c.psi <- ctrl$tuning.psi))

    ## need (r0, r, scale, x, c.psi,c.chi, bb)
    r0 <- obj$init$resid
    r <- resid(obj)
    scale <- obj$scale
    if (is.null(x))  x <- model.matrix(obj)
    bb <- 1/2 ## this is always 1/2 for S estimates by convention
### --- start code from .vcov.MM ---
    ## scaled residuals
    n <- length(r)
    stopifnot(n == length(r0), is.matrix(x), n == nrow(x))
    p <- ncol(x)
    r.s	 <- r / scale # final   scaled residuals
    r0.s <- r0 / scale # initial scaled residuals
    w  <- Mpsi(r.s, cc = c.psi, psi = psi, deriv = 1)
    w0 <- Mchi(r0.s, cc = c.chi, psi = chi, deriv = 1)
    ## FIXME for multivariate y :
    x.wx <- crossprod(x, x * w)
    if(inherits(A <- tryCatch(solve(x.wx) * scale,
			      error=function(e)e), "error")) {
	warning("X'WX is almost singular. Consider rather using cov = \".vcov.w\"")
	A <- tryCatch(solve(x.wx, tol = 0) * scale, error=function(e)e)
	if(inherits(A, "error"))
	    stop("X'WX is singular. Rather use cov = \".vcov.w\"")
    }
    a <- A %*% (crossprod(x, w * r.s) / mean(w0 * r0.s))
    w <- Mpsi( r.s, cc = c.psi, psi = psi)

    ## 3) now the standard part  (w, x, r0.s,  n, A,a, c.chi, bb)
    w0 <- Mchi(r0.s, cc = c.chi, psi = chi)
    Xww <- crossprod(x, w*w0)
    u1 <- A %*% crossprod(x, x * w^2) %*% (n * A)
    u2 <- a %*% crossprod(Xww, A)
    u3 <- A %*% tcrossprod(Xww, a)
    u4 <- mean(w0^2 - bb^2) * tcrossprod(a)

    ## list(cov = matrix((u1 - u2 - u3 + u4)/n, p, p),
    ##      wt = w / r.s, a = a)
### --- end code from .vcov.MM ---
    ret <- (u1 - u2 - u3 + u4)/n

    ## this might not be a positive definite matrix
    ## check eigenvalues (symmetric: ensure non-complex)
    ev <- eigen(ret, symmetric = TRUE)
    if (any(neg.ev <- ev$values < 0)) { ## there's a problem
	posdef.meth <- match.arg(posdef.meth)
	if(ctrl$trace.lev)
	    message("fixing ", sum(neg.ev),
		    " negative eigen([",p,"])values")
	Q <- ev$vectors
	switch(posdef.meth,
	       "orig" = {
		   ## remove negative eigenvalue:
		   ## transform covariance matrix into eigenbasis
		   levinv <- solve(Q)
		   cov.eb <- levinv %*% ret %*% Q
		   ## set vectors corresponding to negative ev to zero
		   cov.eb[, neg.ev] <- 0
		   ## cov.eb[cov.eb < 1e-16] <- 0
		   ## and transform back
		   ret <- Q %*% cov.eb %*% levinv
	       },
	       "posdefify" = {
		   ## Instead of using	require("sfsmisc") and
		   ## ret <- posdefify(ret, "someEVadd",eigen.m = ev,eps.ev = 0)
		   lam <- ev$values
		   lam[neg.ev] <- 0
		   o.diag <- diag(ret)# original one - for rescaling
		   ret <- Q %*% (lam * t(Q)) ## == Q %*% diag(lam) %*% t(Q)
		   ## rescale to the original diagonal values
		   ## D <- sqrt(o.diag/diag(ret))
		   ## where they are >= 0 :
		   D <- sqrt(pmax.int(0, o.diag)/diag(ret))
		   ret[] <- D * ret * rep(D, each = p) ## == diag(D) %*% m %*% diag(D)
	       },
	       stop("invalid 'posdef.meth': ", posdef.meth))
    }
    attr(ret,"weights") <- w / r.s
    attr(ret,"eigen") <- ev
    ret
}## end{.vcov.avar1}

lmrob..M..fit <- function (x=obj$x, y=obj$y, beta.initial=obj$coef,
                           scale=obj$scale, control=obj$control, obj,
                           mf = obj$model
                           ##   ^^^^^^^^^ not model.frame(obj) to avoid errors.
                           )
{
    c.psi <- .psi.conv.cc(control$psi, control$tuning.psi)
    ipsi <- .psi2ipsi(control$psi)
    stopifnot(is.matrix(x))
    n <- nrow(x)
    p <- ncol(x)
    if (is.null(y) && !is.null(obj$model))
        y <- model.response(obj$model, "numeric")
    stopifnot(length(y) == n,
              length(c.psi) > 0, c.psi >= 0,
              scale >= 0, length(beta.initial) == p)

    ret <- .C(R_lmrob_MM,
              x = as.double(x),
              y = as.double(y),
              n = as.integer(n),
              p = as.integer(p),
              beta.initial = as.double(beta.initial),
              scale = as.double(scale),
              coefficients = double(p),
              residuals = double(n),
              iter = as.integer(control$max.it),
              c.psi = as.double(c.psi),
              ipsi = as.integer(ipsi),
              loss = double(1),
              rel.tol = as.double(control$rel.tol),
              converged = logical(1),
              trace.lev = as.integer(control$trace.lev),
              mts = as.integer(control$mts),
              ss =  .convSs(control$subsampling)
              )[c("coefficients",  "scale", "residuals", "loss", "converged", "iter")]
    ## FIXME?: Should rather warn *here* in case of non-convergence
    ret$fitted.values <- drop(x %*% ret$coefficients)
    names(ret$coefficients) <- colnames(x)
    names(ret$residuals) <- rownames(x)
    ret$rweights <- lmrob.rweights(ret$residuals, scale, control$tuning.psi, control$psi)
    if (!grepl('M$', control$method)) {
        ## update control$method if it's not there already
        control$method <- paste0(control$method, 'M')
    }
    ret$control <- control
    if (!missing(obj)) {
        if (!is.null(obj$call)) {
            ret$call <- obj$call
            ret$call$method <- control$method
        }
        if (control$method %in% c('SM', 'MM')) {
            ret$init.S <- obj
        } else {
            ret$init <-
                obj[names(obj)[na.omit(match(
                    c("coefficients","scale", "residuals", "loss", "converged",
                      "iter", "rweights", "fitted.values", "control", "ostats",
                      "init.S", "init", "kappa", "tau"), names(obj)))]]
            class(ret$init) <- 'lmrob'
            ret <- c(ret,
                     obj[names(obj)[na.omit(match(
                         c("df.residual", "degree.freedom",
                           "xlevels", "terms", "model", "x", "y",
                           "na.action", "contrasts", "MD"), names(obj)))]])
        }
        ret$qr <- qr(x * sqrt(ret$rweights))
        ret$rank <- ret$qr$rank
        ## if there is a covariance matrix estimate available in obj
        ## update it, if possible, else replace it by the default .vcov.w
	if (!is.null(obj$cov)) {
	    if (!control$method %in% c('SM', 'MM') &&
		ret$control$cov == '.vcov.avar1')
		ret$control$cov <- '.vcov.w'
	    lf.cov <- if (!is.function(ret$control$cov))
		get(ret$control$cov, mode='function') else ret$control$cov
	    ret$cov <- lf.cov(ret, x)
	}
        if (!is.null(obj$assign)) ret$assign <- obj$assign
    }
    class(ret) <- "lmrob"
    if (control$method %in% control$compute.outlier.stats & !missing(obj))
        ret$ostats <- outlierStats(ret, x, control)
    ret
}


lmrob.S <- function (x, y, control, trace.lev = control$trace.lev, mf = NULL)
{
    if (!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    nResample <- as.integer(control$nResample)
    groups <- as.integer(control$groups)
    nGr <- as.integer(control$n.group)
    large_n <- (n > control$fast.s.large.n)
    if (large_n) {
        if (nGr <= p)
            stop("'control$n.group' must be larger than 'p' for 'large_n' algorithm")
        if (nGr * groups > n)
            stop("'groups * n.group' must be smaller than 'n' for 'large_n' algorithm")
        if (nGr <= p + 10) ## FIXME (be smarter ..)
            warning("'control$n.group' is not much larger than 'p', probably too small")
    }
    if (length(seed <- control$seed) > 0) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            seed.keep <- get(".Random.seed", envir = .GlobalEnv,
                             inherits = FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
        }
        assign(".Random.seed", seed, envir = .GlobalEnv) ## why not set.seed(seed)
    }

    bb <- as.double(control$bb)
    c.chi <- .psi.conv.cc(control$psi, control$tuning.chi)
    best.r <- as.integer(control$best.r.s)
    stopifnot(length(c.chi) > 0, c.chi >= 0, length(bb) > 0,
              length(best.r) > 0, best.r >= 1, length(y) == n, n > 0)

    b <- .C(R_lmrob_S,
            x = as.double(x),
            y = as.double(y),
            n = as.integer(n),
            p = as.integer(p),
            nResample = nResample,
            scale = double(1),
            coefficients = double(p),
            as.double(c.chi),
            .psi2ipsi(control$psi),
            bb,
            best_r = best.r,
            groups = groups,
            n.group = nGr,
            k.fast.s = as.integer(control$k.fast.s),
            k.iter = as.integer(control$k.max),
            maxit.scale = as.integer(control$maxit.scale),
            refine.tol = as.double(control$refine.tol),
            inv.tol = as.double(control$solve.tol),
            converged = logical(1),
            trace.lev = as.integer(trace.lev),
            mts = as.integer(control$mts),
            ss = .convSs(control$subsampling),
            fast.s.large.n = as.integer(if (large_n) control$fast.s.large.n else n+1)
            ## avoids the use of NAOK = TRUE for control$fast.s.large.n == Inf
            )[c("coefficients", "scale", "k.iter", "converged")]
    scale <- b$scale
    if (scale < 0)
	stop("C function R_lmrob_S() exited prematurely")
    if (scale == 0)
	warning("S-estimated scale == 0:  Probably exact fit; check your data")
    ## FIXME: get 'res'iduals from C

    b$fitted.values <- x %*% b$coefficients
    b$residuals <- drop(y - b$fitted.values)
    names(b$coefficients) <- colnames(x)
    names(b$residuals) <- rownames(x)
    ## robustness weights
    b$rweights <- lmrob.rweights(b$residuals, scale, control$tuning.chi, control$psi)
    ## set method argument in control
    control$method <- 'S'
    b$control <- control
    ## add call if called from toplevel
    if (identical(parent.frame(), .GlobalEnv))
        b$call <- match.call()
    class(b) <- 'lmrob.S'
    if ("S" %in% control$compute.outlier.stats)
        b$ostats <- outlierStats(b, x, control)
    b
}

lmrob..D..fit <- function(obj, x=obj$x, control = obj$control,
                          mf = obj$model)
{
    if (is.null(control)) stop('lmrob..D..fit: control is missing')
    if (!obj$converged)
        stop('lmrob..D..fit: prior estimator did not converge, stopping')
    if (is.null(x)) x <- model.matrix(obj)
    w <- obj$rweights
    if (is.null(w)) stop('lmrob..D..fit: robustness weights undefined')
    if (is.null(obj$residuals)) stop('lmrob..D..fit: residuals undefined')
    r <- obj$residuals
    psi <- control$psi
    if (is.null(psi)) stop('lmrob..D..fit: parameter psi is not defined')
    c.psi <- .psi.conv.cc(psi, if (control$method %in% c('S', 'SD'))
                           control$tuning.chi else control$tuning.psi)
    if (!is.numeric(c.psi)) stop('lmrob..D..fit: parameter tuning.psi is not numeric')

    obj$init <- obj[names(obj)[na.omit(match(
        c("coefficients","scale", "residuals", "loss", "converged", "iter", "ostats",
	  "rweights", "fitted.values", "control", "init.S", "init"), names(obj)))]]
    obj$init.S <- NULL

    if (is.null(obj$kappa))
        obj$kappa <- lmrob.kappa(obj, control)
    kappa <- obj$kappa
    if (is.null(obj$tau))
        obj$tau <- lmrob.tau(obj, x, control)
    tau <- obj$tau

    ## get starting value for root search (to keep breakdown point!)
    scale.1 <- sqrt(sum(w * r^2) / kappa / sum(tau^2*w))
    ret <- .C(R_find_D_scale,
              r = as.double(r),
              kappa = as.double(kappa),
              tau = as.double(tau),
              length = as.integer(length(r)),
              scale = as.double(scale.1),
              c = as.double(c.psi),
              ipsi = .psi2ipsi(psi),
              type = 3L, ## dt1 as only remaining option
              rel.tol = as.double(control$rel.tol),
              k.max = as.integer(control$k.max),
              converged = logical(1))[c("converged", "scale")]
    obj$scale <- if(ret$converged) ret$scale else NA
    obj$converged <- ret$converged

    if (!grepl('D$', control$method)) {
        ## append "D"  to control$method if it's not there already
        method <- control$method
        if (method == 'MM') method <- 'SM'
        control$method <- paste0(method, 'D')
    }
    ## update call
    if (!is.null(obj$call)) obj$call$method <- control$method
    obj$control <- control
    class(obj) <- "lmrob"

    ## if there is a covariance matrix estimate available in obj
    ## update it, if possible, else replace it by the default
    ## .vcov.w
    if (!is.null(obj$cov)) {
        if (control$cov == '.vcov.avar1')
            control$cov <- '.vcov.w'

        lf.cov <- if (!is.function(control$cov))
            get(control$cov, mode='function') else control$cov
        obj$cov <- lf.cov(obj, x)
    }
    if (control$method %in% control$compute.outlier.stats)
        obj$ostats <- outlierStats(obj, x, control)

    obj
}

globalVariables(c("psi", "wgt", "r"), add=TRUE) ## <- lmrob.E( <expr> )

lmrob.kappa <- function(obj, control = obj$control)
{
    if (is.null(control)) stop('control is missing')
    if (control$method %in% c('S', 'SD')) control$tuning.psi <- control$tuning.chi

    fun.min <- function(kappa) lmrob.E(psi(r)*r - kappa*wgt(r), control = control)
    uniroot(fun.min, c(0.1, 1))$root
}

## "FIXME" How to get \hat{tau} for a simple *M* estimate here ??
## lmrob.tau() is used in lmrob..D..fit()
lmrob.tau <- function(obj, x=obj$x, control = obj$control, h, fast = TRUE)
{
    if(is.null(control)) stop("'control' is missing")
    if(missing(h))
	h <- if (is.null(obj$qr))
	    lmrob.leverages(x, obj$rweights)
	else
	    lmrob.leverages(wqr = obj$qr)

    ## speed up: use approximation if possible
    if (fast && !control$method %in% c('S', 'SD')) {
        c.psi <- control$tuning.psi
        tfact <- tcorr <- NA
        switch(control$psi,
               optimal = if (isTRUE(all.equal(c.psi, 1.060158))) {
                   tfact <- 0.94735878
                   tcorr <- -0.09444537
               },
               bisquare = if (isTRUE(all.equal(c.psi, 4.685061))) {
                   tfact <- 0.9473684
                   tcorr <- -0.0900833
               },
               welsh = if (isTRUE(all.equal(c.psi, 2.11))) {
                   tfact <- 0.94732953
                   tcorr <- -0.07569506
               },
               ggw = if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.95, NA)))) {
                   tfact <- 0.9473787
                    tcorr <- -0.1143846
               } else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.95, NA)))) {
                   tfact <- 0.94741036
                   tcorr <- -0.08424648
               },
               lqq = if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.95, NA)))) {
                   tfact <- 0.94736359
                   tcorr <- -0.08594805
               },
               hampel = if (isTRUE(all.equal(c.psi, c(1.35241275, 3.15562975, 7.212868)))) {
                   tfact <- 0.94739770
                   tcorr <- -0.04103958
               },
           {})
        if (!is.na(tfact))
            return(sqrt(1 - tfact*h) * (tcorr*h + 1))
    }
    ## else "non-fast" -- need to compute the integrals :

    ## kappa
    kappa <- if(is.null(obj$kappa)) lmrob.kappa(obj, control) else obj$kappa
    ## local variables
    ## n <- length(h)
    ## set psi and cpsi
    psi <- control$psi
    if (is.null(psi)) stop('parameter psi is not defined')
    cpsi <- if (control$method %in% c('S', 'SD'))
        control$tuning.chi else control$tuning.psi
    cpsi <- .psi.conv.cc(psi, cpsi)# has its test
    ipsi <- .psi2ipsi(psi)

    ## constant for stderr of u_{-i} part and other constants
    inta <- function(r) .Mpsi(r, cpsi, ipsi)^2          * dnorm(r)
    intb <- function(r) .Mpsi(r, cpsi, ipsi, deriv = 1) * dnorm(r)
    ## intc <- function(r) .Mpsi(r, cpsi, ipsi) * r        * dnorm(r)
                                        # changed from psi/e to psi*e
    ta <- integrate(inta, -Inf,Inf)$value
    tb <- integrate(intb, -Inf,Inf)$value
    ## tE <- integrate(intc, -Inf,Inf)$value

    ## calculate tau for unique h
    hu <- unique(h)
    nu <- length(hu)

    ## Initialize tau vector
    tau <- numeric(length=nu)

    tc <- ta/tb^2
    ## --- Gauss-Hermite integration
    gh <- ghq(control$numpoints)
    ghz <- gh$nodes
    ghw <- gh$weights
    ## Calulate each tau_i
    for (i in 1:nu) {
        ## stderr of u_{-i} part
        s <- sqrt(tc*(hu[i]-hu[i]^2))
        tc2 <- hu[i]/tb
        ## function to be integrated
        fun <- function(w, v, sigma.i) {
	    t <- (v - tc2*.Mpsi(v, cpsi, ipsi) + w*s)/sigma.i
	    psi.t <- .Mpsi(t, cpsi, ipsi)
	    (psi.t*t - kappa*psi.t/t) * dnorm(v)*dnorm(w)
        }
        ## integrate over w
        wint <- function(v, sigma.i) {
            ## sapply(v,function(v.j) integrate(fun,-Inf,Inf,v.j,sigma.i)$value)
            sapply(v, function(v.j) sum(fun(ghz, v.j, sigma.i)*ghw))
        }
        ## integrate over v
        vint <- function(sigma.i) {
            ## integrate(wint,-Inf,Inf,sigma.i)$value
            sum(wint(ghz, sigma.i)*ghw)
        }

        ## find tau
        tau[i] <- uniroot(vint, c(if (hu[i] < 0.9) 3/20 else 1/16, 1.1))$root
    }

    tau[match(h, hu)]
}

lmrob.tau.fast.coefs <- function(cc, psi) {
    ## function that calculates the coefficients for 'fast' mode of lmrob.tau
    ctrl <- lmrob.control(tuning.psi = cc, psi = psi)
    levs <- seq(0, 0.8, length.out = 80)
    ## calculate taus
    taus <- lmrob.tau(list(), control=ctrl, h=levs, fast=FALSE)
    ## calculate asymptotic approximation of taus
    ta <- lmrob.E(psi(r)^2,  ctrl, use.integrate = TRUE)
    tb <- lmrob.E(psi(r, 1), ctrl, use.integrate = TRUE)
    tfact <- 2 - ta/tb^2
    taus.0 <- sqrt(1 - tfact * levs)
    ## calculate correction factor
    tcorr <- coef(lmrob(taus / taus.0 - 1 ~ levs - 1))
    c(tfact = tfact, tcorr = tcorr)
}

lmrob.hatmatrix <- function(x, w = rep(1, NROW(x)), wqr = qr(sqrt(w) * x))
{
    tcrossprod(qr.Q(wqr))
}

lmrob.leverages <- function(x, w = rep(1, NROW(x)), wqr = qr(sqrt(w) * x))
{
    if (missing(wqr) && !is.matrix(x)) x <- as.matrix(x)
    ## Faster than computing the whole hat matrix, and then diag(.) :
    ## == diag(lmrob.hatmatrix(x, w, ...))
    pmin(1, rowSums(qr.Q(wqr)^2))
}

##' psi |--> ipsi \in \{0,1,...6} : integer codes used in C
.psi2ipsi <- function(psi)
{
    psi <- .regularize.Mpsi(psi, redescending=FALSE)
    i <- match(psi, c(
	'huber', 'bisquare', 'welsh', 'optimal',
	## 0	    1	        2	 3
	'hampel', 'ggw', 'lqq'
	## 4	    5	   6
	))
    if(is.na(i)) stop("internal logic error in psi() function name: ", psi,
		      "  Please report!")
    i - 1L
}

##' Given psi() fn (as string), possibly convert the tuning-constant vector cc
##' such that it "fits" to psi()
.psi.conv.cc <- function(psi, cc)
{
    if (!is.character(psi) || length(psi) != 1)
        stop("argument 'psi' must be a string (denoting a psi function)")
    if(!is.numeric(cc))
        stop("tuning constant 'cc' is not numeric")

    ## "FIXME": For (ggw, lqq) this is much related to  .psi.const() below
    switch(tolower(psi),
           'ggw' = {
               ## Input: 4 parameters, (minimal slope, b, efficiency, breakdown point)
               ## Output 'k': either k in {1:6} or  k = c(0, k[2:5])
               if (isTRUE(all.equal(cc, c(-.5, 1, 0.95, NA)))) return(1)
               else if (isTRUE(all.equal(cc, c(-.5, 1, 0.85, NA)))) return(2)
               else if (isTRUE(all.equal(cc, c(-.5, 1.0, NA, 0.5)))) return(3)
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA)))) return(4)
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, 0.85, NA)))) return(5)
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5)))) return(6)
               else if (length(cc) == 5 && cc[1] == 0 ||
                        (length(cc <- attr(cc, 'constants')) == 5 && cc[1] == 0))
                   return(cc)
               else stop('Coefficients for ',psi,' function incorrectly specified.\n',
                         'Use c({0, } minimal slope, b, efficiency, breakdown point)')
           },
           'lqq' = {
               ## Input: 4 parameters, (minimal slope, b, efficiency, breakdown point)
               ## Output: k[1:3] = (b, c, s)
               if (isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA))))
                   return(c(1.4734061, 0.9822707, 1.5))
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5))))
                   return(c(0.4015457, 0.2676971, 1.5))
               else if (length(cc) == 3 || length(cc <- attr(cc, 'constants')) == 3)
                   return(cc)
               else stop('Coefficients for ',psi,' function incorrectly specified.\n',
                         'Use c(minimal slope, b, efficiency, breakdown point) or k[1:3]')
           },
           'hampel' = {
               ## just check length of coefficients
               if (length(cc) != 3)
                   stop('Coef. for Hampel psi function not of length 3')
           }, {
               ## otherwise: should have length 1
               if (length(cc) != 1)
                   stop('Coef. for psi function ', psi,' not of length 1')
           })

    return(cc)
}
##' @title For GGW's psi(), find x with minimal slope, and the min.slope
##' @param a "scale" of GGW's psi
##' @param b exponent of GGW's psi
##' @param c "huber-cutoff" of GGW's psi
##' @param ... further arguments passed to optimize()
##' @return the return value of optimize():  list(minimum, objective)
##' @author Manuel Kohler and Martin Maechler
lmrob.ggw.mxs <- function(a, b, c, ...) {
    ipsi <- .psi2ipsi('ggw')
    ccc <- c(0, a, b, c, 1) ## == .psi.conv.cc('ggw', cc=c(0, a, b, c, 1))
    optimize(.Mpsi, c(c, max(a+b+2*c, 0.5)), ccc=ccc, ipsi=ipsi,
	     deriv = 1, ...)
}

lmrob.ggw.ms <- function(a, b, c, ...) ## find minimal slope
    lmrob.ggw.mxs(a, b, c, ...)[["objective"]]

lmrob.ggw.finda <- function(ms, b, c, ...) ## find a constant
{
    val <- uniroot(function(a) lmrob.ggw.ms(1/a, b, c) - ms,
                   c(200, if (b > 1.4) 1/400 else if (b > 1.3) 1/50 else 1/20),
                   ...)
    1/val$root
}

lmrob.ggw.ac <- function(a, b, c) ## calculate asymptotic efficiency
{
    ipsi <- .psi2ipsi('ggw')
    ccc <- c(0, a, b, c, 1)
    lmrob.E(.Mpsi(r, ccc, ipsi, deriv=1), use.integrate = TRUE)^2 /
    lmrob.E(.Mpsi(r, ccc, ipsi) ^2,       use.integrate = TRUE)
}

lmrob.ggw.bp <- function(a, b, c, ...) { ## calculate kappa
    ipsi <- .psi2ipsi('ggw')
    abc <- c(0, a, b, c)
    nc <- integrate(.Mpsi, 0, Inf, ccc = c(abc, 1), ipsi=ipsi, ...)$value
    lmrob.E(.Mchi(r, ccc = c(abc, nc), ipsi), use.integrate = TRUE)
}

.psi.ggw.findc <- function(ms, b, eff = NA, bp = NA) {
    ## find c by eff for bp
    c <- if (!is.na(eff)) {
        if (!is.na(bp))
            warning('tuning constants for ggw psi: both eff and bp specified, ignoring bp')
        ## find c by eff
        uniroot(function(x) lmrob.ggw.ac(lmrob.ggw.finda(ms, b, x), b, x) - eff,
                c(0.15, if (b > 1.61) 1.4 else 1.9))$root
    } else {
        if (is.na(bp))
	    stop("neither breakdown point 'bp' nor efficiency 'eff' specified")
        ## find c by bp
        uniroot(function(x) lmrob.ggw.bp(lmrob.ggw.finda(ms, b, x), b, x) - bp,
                c(0.08, if (ms < -0.4) 0.6 else 0.4))$root
    }

    a <- lmrob.ggw.finda(ms, b, c)
    ipsi <- .psi2ipsi('ggw')
    ccc <- c(0, a, b, c, 1)
    nc <- integrate(.Mpsi, 0, Inf, ccc=ccc, ipsi=ipsi)$value
    ## return
    c(0, a, b, c, nc)
}

lmrob.efficiency <-  function(psi, cc, ...) {
    ipsi <- .psi2ipsi(psi)
    ccc <- .psi.conv.cc(psi, cc=cc)

    integrate(function(x) .Mpsi(x, ccc=ccc, ipsi=ipsi, deriv=1)*dnorm(x),
	      -Inf, Inf, ...)$value^2 /
    integrate(function(x) .Mpsi(x, ccc=ccc, ipsi=ipsi)^2 *dnorm(x),
	      -Inf, Inf, ...)$value
}

lmrob.bp <- function(psi, cc, ...)
  integrate(function(x) Mchi(x, cc, psi)*dnorm(x), -Inf, Inf, ...)$value

##' @title Find tuning constant for  "lqq"  psi function
##' @param cc numeric vector =  c(min_slope, b, eff, bp);
##'      typically 'eff' or 'bp' are NA and will be computed
##' @param interval for finding 'c' via uniroot()
##' @param subdivisions for integrate()
##' @param rel.tol relative and
##' @param abs.tol absolute tolerance for integrate()
##' @param tol relative tolerance for uniroot()
##' @param maxiter maximal number of iterations for uniroot()
##' @return constants for c function: c(b*c, c, s = 1 - min_slope)
##' @author Manuel Koller and Martin Maechler
.psi.lqq.findc <-
    function(cc, interval = c(0.1, 4),
             subdivisions = 100L,
             rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
             tol = .Machine$double.eps^0.25, maxiter = 1000)
{
    t.fun <- if (!is.na(cc[3])) {
	if (!is.na(cc[4]))
	    warning('tuning constants for lqq psi: both eff and bp specified, ignoring bp')
	## find c by b, s and eff
	function(c)
	    lmrob.efficiency('lqq', c(cc[2]*c, c, 1-cc[1]),
                             subdivisions=subdivisions,
                             rel.tol=rel.tol, abs.tol=abs.tol) - cc[3]
    } else {
	if (is.na(cc[4]))
	    stop('Error: neither breakdown point nor efficiency specified')
	function(c)
	    lmrob.bp('lqq', c(cc[2]*c, c, 1-cc[1]),
                     subdivisions=subdivisions,
                     rel.tol=rel.tol, abs.tol=abs.tol) - cc[4]
    }
    c. <- tryCatch(uniroot(t.fun, interval=interval, tol=tol, maxiter=maxiter)$root,
		   error=function(e)e)
    if (inherits(c., 'error'))
        stop('.psi.lqq.findc: unable to find constants for psi function')
    else c(cc[2]*c., c., 1-cc[1])
}

##' For ("ggw", "lqq"), if  cc is not one of the predefined ones,
##'
##' compute the tuning constants numerically, from the given specs (eff / bp):
.psi.const <- function(cc, psi)
{
    switch(psi,
           "ggw" = { ## only calculate for non-standard coefficients
               if (!(isTRUE(all.equal(cc, c(-.5, 1,   0.95, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1,   0.85, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1,   NA,  0.5))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, 0.85, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5))))) {
		   attr(cc, 'constants') <-
			.psi.ggw.findc(ms=cc[1], b=cc[2], eff=cc[3], bp=cc[4])
               }
           },
           "lqq" = { ## only calculate for non-standard coefficients
               if (!(isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5))))) {
                   attr(cc, 'constants') <- .psi.lqq.findc(cc)
               }
           },
           stop("method for psi function ", psi, " not implemented"))
    cc
}

Mpsi <- function(x, cc, psi, deriv=0) {
    x[] <- .Call(R_psifun, x, .psi.conv.cc(psi, cc), .psi2ipsi(psi), deriv)
    x
}
.Mpsi <- function(x, ccc, ipsi, deriv=0) .Call(R_psifun, x, ccc, ipsi, deriv)

Mchi <- function(x, cc, psi, deriv=0) {
    x[] <- .Call(R_chifun, x, .psi.conv.cc(psi, cc), .psi2ipsi(psi), deriv)
    x
}
.Mchi <- function(x, ccc, ipsi, deriv=0) .Call(R_chifun, x, ccc, ipsi, deriv)

Mwgt <- function(x, cc, psi) {
    x[] <- .Call(R_wgtfun, x, .psi.conv.cc(psi, cc), .psi2ipsi(psi))
    x
}
.Mwgt <- function(x, ccc, ipsi) .Call(R_wgtfun, x, ccc, ipsi)

## only for nlrob() -- and to use instead of MASS:::psi.huber etc:
## returns a *function* a la  psi.huber() :
.Mwgt.psi1 <- function(psi, cc = .Mpsi.tuning.default(psi)) {
    ipsi <- .psi2ipsi(psi)
    ccc <- .psi.conv.cc(psi, cc)
    ## return function *closure* :
    function(x, deriv = 0)
    if(deriv) .Mpsi(x, ccc, ipsi, deriv=deriv) else .Mwgt(x, ccc, ipsi)
}

##' The normalizing constant for  rho(.) <--> rho~(.)
MrhoInf <- function(cc, psi) {
    cc <- .psi.conv.cc(psi, cc)
    .Call(R_rho_inf, cc, .psi2ipsi(psi))
}
.MrhoInf <- function(ccc, ipsi) .Call(R_rho_inf, ccc, ipsi)


lmrob.rweights <- function(resid, scale, cc, psi, eps = 16 * .Machine$double.eps) {
    if (scale == 0) { ## exact fit
	m <- max(ar <- abs(resid))
	if(m == 0) numeric(seq_len(ar)) else
	as.numeric(ar < eps * m)# 1 iff res ~= 0
    } else
	Mwgt(resid / scale, cc, psi)
}

lmrob.E <- function(expr, control, dfun = dnorm, use.integrate = FALSE, obj, ...)
{
    expr <- substitute(expr)
    if (missing(control) && !missing(obj))
        control <- obj$control

    lenvir <-
      if (!missing(control)) {
        psi <- control$psi
        if (is.null(psi)) stop('parameter psi is not defined')
	c.psi <- control[[if (control$method %in% c('S', 'SD'))
			  "tuning.chi" else "tuning.psi"]]
        if (!is.numeric(c.psi)) stop('tuning parameter (chi/psi) is not numeric')

        list(psi = function(r, deriv = 0) Mpsi(r, c.psi, psi, deriv),
             chi = function(r, deriv = 0) Mchi(r, c.psi, psi, deriv), ## change?
             wgt = function(r) Mwgt(r, c.psi, psi)) ## change?

      } else list()

    pf <- parent.frame()
    FF <- function(r)
	eval(expr, envir = c(list(r = r), lenvir), enclos = pf) * dfun(r)
    if (isTRUE(use.integrate)) {
	integrate(FF, -Inf,Inf, ...)$value
    ## This would be a bit more accurate .. *AND* faster notably for larger 'numpoints':
    ## } else if(use.integrate == "GQr") {
    ##     require("Gqr")# from R-forge [part of lme4 project]
    ##     ## initialize Gauss-Hermite Integration
    ##     GH <- GaussQuad(if(is.null(control$numpoints)) 13 else control$numpoints,
    ##                     "Hermite")
    ##     ## integrate
    ##     F. <- function(r) eval(expr, envir = c(list(r = r), lenvir), enclos = pf)
    ##     sum(GH$weights * F.(GH$knots))
    } else {
	## initialize Gauss-Hermite Integration
	gh <- ghq(if(is.null(control$numpoints)) 13 else control$numpoints)
	## integrate
	sum(gh$weights * FF(gh$nodes))
    }
}

ghq <- function(n = 1, modify = TRUE) {
    ## Adapted from gauss.quad in statmod package
    ## which itself has been adapted from Netlib routine gaussq.f
    ## Gordon Smyth, Walter and Eliza Hall Institute

    n <- as.integer(n)
    if(n<0) stop("need non-negative number of nodes")
    if(n==0) return(list(nodes=numeric(0), weights=numeric(0)))
    ## i <- seq_len(n) # 1 .. n
    i1 <- seq_len(n-1L)

    muzero <- sqrt(pi)
    ## a <- numeric(n)
    b <- sqrt(i1/2)

    A <- numeric(n*n)
    ## A[(n+1)*(i-1)+1] <- a # already 0
    A[(n+1)*(i1-1)+2] <- b
    A[(n+1)*i1] <- b
    dim(A) <- c(n,n)
    vd <- eigen(A,symmetric=TRUE)
    n..1 <- n:1L
    w <- vd$vectors[1, n..1]
    w <- muzero * w^2
    x <- vd$values[n..1] # = rev(..)
    list(nodes=x, weights= if (modify) w*exp(x^2) else w)
}

.convSs <- function(ss)
    switch(ss,
           "simple"= 0L,
           "nonsingular"= 1L,
           stop("unknown setting for parameter ss"))


outlierStats <- function(object, x = object$x,
                         control = object$control,
                         epsw = control$eps.outlier,
                         epsx = control$eps.x,
                         warn.limit.reject = control$warn.limit.reject,
                         warn.limit.meanrw = control$warn.limit.meanrw
                         ) {
    ## look at all the factors in the model and count
    ## for each level how many observations were rejected.
    ## Issue a warning if there is any level where more than
    ## warn.limit.reject observations were rejected or
    ## the mean robustness weights was <= warn.limit.meanrw

    rw <- object$rweights
    ##    ^^^^^^^^^^^^^^^ not weights(..., type="robustness") as we
    ##                    don't want naresid() padding here.
    if (is.function(epsw)) epsw <- epsw(nobs(object))
    if (!is.numeric(epsw) || length(epsw) != 1)
        stop("'epsw' must be numeric(1) or a function of nobs(obj.) which returns a numeric(1)")
    rj <- abs(rw) < epsw
    if (NROW(x) != length(rw))
        stop("number of rows in 'x' and length of 'object$rweights' must be the same")
    if (is.function(epsx)) epsx <- epsx(max(abs(x)))
    if (!is.numeric(epsx) || length(epsx) != 1)
        stop("'epsx' must be numeric(1) or a function of max(abs(x)) which returns a numeric(1)")
    xnz <- abs(x) > epsx

    cc <- function(idx) {
        nnz <- sum(idx) ## <- if this is zero, 'Ratio' and 'Mean.RobWeight' will be NaN
        Fr <- sum(rj[idx])
	c(N.nonzero = nnz,
	  N.rejected = Fr,
	  Ratio = Fr / nnz,
	  Mean.RobWeight = mean(rw[idx]))
    }

    report <- t(apply(cbind(Overall=TRUE, xnz[, colSums(xnz) < NROW(xnz)]), 2, cc))

    shout <- FALSE # should we "shout"? -- scalar logical, never NA
    lbr <- rep.int(FALSE, nrow(report))
    if (!is.null(warn.limit.reject)) {
	lbr <- report[, "Ratio"] >= warn.limit.reject
	shout <- any(lbr & !is.na(lbr))
    }
    if (!is.null(warn.limit.meanrw)) {
	lbr <- lbr | report[, "Mean.RobWeight"] <= warn.limit.meanrw
	shout <- shout || any(lbr & !is.na(lbr))
    }
    if (shout) {
        nbr <- rownames(report)[lbr]
        attr(report, "warning") <- paste("Possible local breakdown of",
                                         paste0("'", nbr, "'", collapse=", "))
        warning("Detected possible local breakdown of ", control$method, "-estimate in ",
                if (length(nbr) > 1) paste(length(nbr), "coefficients") else "coefficient",
                " ", paste0("'", nbr, "'", collapse=", "), ".",
                if ("KS2014" %in% control$setting) "" else
                "\nUse lmrob argument 'setting=\"KS2014\"' to avoid this problem."
                )
    }

    report
}
