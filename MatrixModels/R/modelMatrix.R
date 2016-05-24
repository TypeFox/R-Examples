####---- This was part of ../../Matrix/R/spModels.R -- till 2010-07-25

model.Matrix <- function(object, data = environment(object),
		 contrasts.arg = NULL, xlev = NULL,
			 sparse = FALSE, drop.unused.levels = FALSE, ...)
{
    if(sparse) {
	m <- sparse.model.matrix(object, data=data, contrasts.arg=contrasts.arg,
				 drop.unused.levels=drop.unused.levels, xlev=xlev,
				 ...)
	new("dsparseModelMatrix",  m, ## dropping attributes ?
	    assign = attr(m, "assign"),
	    contrasts = if(is.null(ctr <- attr(m,"contrasts")))list() else ctr)
    } else {
	## as standard	model.matrix() but producing a	"ddenseModelMatrix":
	m <- model.matrix(object, data=data,
			  contrasts.arg=contrasts.arg, xlev=xlev, ...)
	new("ddenseModelMatrix", as(m, "dgeMatrix"),
	    assign = attr(m, "assign"),
	    contrasts = if(is.null(ctr <- attr(m,"contrasts")))list() else ctr)
    }
}



### Keep this namespace-hidden: Would need to return a classed object

## FIXME: still test this function for both methods, since currently
## ----- both  dgCMatrix_cholsol and  dgCMatrix_qrsol are only called from here!
lm.fit.sparse <- function(x, y, w = NULL, offset = NULL,
			  method = c("qr", "cholesky"),
			  tol = 1e-7, singular.ok = TRUE, order = NULL,
			  transpose = FALSE)
### Fit a linear model, __ given __ a sparse model matrix 'x'
### using a sparse QR or a sparse Cholesky factorization
{
    cld <- getClass(class(x))
    stopifnot(extends(cld, "dsparseMatrix"), is.numeric(y))
## or	  if(!is(x, "dsparseMatrix")) x <- as(x, "dsparseMatrix")
    if(transpose) { tx <- x ; x <- t(x) }
    n <- nrow(x)
    if(NROW(y) != n) stop("incompatible dimensions of (x,y)")
    ny <- NCOL(y)
    if (!is.null(offset)) {
	stopifnot(length(offset) == n)
	y <- y - as.numeric(offset)
    }
    if(ny != 1L) ## FIXME: should not be too much work!
	stop("multivariate, i.e., matrix 'y' is not yet implemented")
    if ((has.w <- !is.null(w))) {
	if(any(w < 0 | is.na(w)))
	    stop("missing or negative weights not allowed")
	if(length(w) != n)
	    stop("weights vector 'w' is of wrong length")

	zero.weights <- any(wis0 <- w == 0)
	if (zero.weights) {
	    save.r <- y
	    save.f <- y
	    save.w <- w
	    ok <- !wis0 # == w != 0
	    i0 <- which(wis0)
	    ok <- which(ok) # (faster when indexing repeatedly)
	    w <- w[ok]
	    x0 <- x[i0, , drop = FALSE]
	    x  <- x[ok, , drop = FALSE]
	    n <- nrow(x)
	    y0 <- if (ny > 1L) y[i0, , drop = FALSE] else y[i0]
	    y  <- if (ny > 1L) y[ok, , drop = FALSE] else y[ok]
	}
	wts <- sqrt(w)
	## keep the unweighted (x,y):
	y. <- y ## x. <- x
	x <- x * wts
	y <- y * wts
    }

    method <- match.arg(method)
    order <- {
	if(is.null(order)) ## recommended default depends on method :
	    if(method == "qr") 3L else 1L
	else as.integer(order) }

    switch(method,
	   "cholesky" = {
	       r <- .solve.dgC.chol(as(if(transpose) tx else t(x), "CsparseMatrix"), y)
	       coef <- r[["coef"]]
	   },
	   "qr" = {
	       coef <-
                   .solve.dgC.qr(if(cld@className %in% c("dtCMatrix", "dgCMatrix")) x
                                 else as(x, "dgCMatrix"),
                                 y, order)
	       ## for now -- FIXME --
	       return(coef)
	   },
	   ## otherwise:
	   stop("unknown method ", dQuote(method))
	   )

    ## FIXME: add names to coef as in lm.wfit(),
    ##		~/R/D/r-devel/R/src/library/stats/R/lm.R
    resid <- if(has.w) r[["resid"]] / wts else r[["resid"]]
    z <- list(coef = coef, weights = w,
	      residuals = resid, fitted.values = y - resid)
    if(has.w && zero.weights) {
	coef[is.na(coef)] <- 0
	f0 <- x0 %*% coef
	if (ny > 1) {
	    save.r[ok, ] <- resid
	    save.r[i0, ] <- y0 - f0
	    save.f[ok, ] <- z$fitted.values
	    save.f[i0, ] <- f0
	}
	else {
	    save.r[ok] <- resid
	    save.r[i0] <- y0 - f0
	    save.f[ok] <- z$fitted.values
	    save.f[i0] <- f0
	}
	z$residuals <- save.r
	z$fitted.values <- save.f
	z$weights <- save.w
    }
    if(!is.null(offset))
	z$fitted.values <- z$fitted.values + offset

    z
}

## allow extra args to be passed to print, notably those
## to printSpMatrix()  [ ../sparseMatrix.R ] :
printModelMat <- function(x, ...)
{
    ## workaround because	 callNextMethod() fails here:
    cat(sprintf("\"%s\": ", class(x)[1]))
    ## (an "intermediate" class) - why exactly? -- callNextMethod()
    print(as(x, "generalMatrix"), ...)
    ## end{workaround}
    p <- length(ass <- x@assign)
    c.ass <- encodeString(ass)
    if(sum(nchar(c.ass))+ p-1 < getOption("width") - 10) ## short enough
	cat("@ assign: ", c.ass,"\n")
    else {
	cat("@ assign:\n"); print(ass)
    }
    cat("@ contrasts:\n"); print(x@contrasts)
    invisible(x)
}

setMethod("print", "modelMatrix", printModelMat)
setMethod("show", "modelMatrix", function(object) printModelMat(object))


setAs("ddenseModelMatrix", "predModule",
      function(from)
  {
      p <- ncol(from)
      new("dPredModule", coef = numeric(p), Vtr = numeric(p),
          X = from, fac = chol(crossprod(from)))
  })

setAs("dsparseModelMatrix", "predModule",
      function(from)
  {
      p <- ncol(from)
      new("sPredModule", coef = numeric(p), Vtr = numeric(p),
          X = from, fac = Cholesky(crossprod(from), LDL = FALSE))
  })

##' Create an respModule, which could be from a derived class such as
##' glmRespMod or nlsRespMod.
##' @title Create a respModule object
##' @param a model frame
##' @param family the optional glm family (glmRespMod only)
##' @param nlenv the nonlinear model evaluation environment (nlsRespMod only)
##' @param nlmod the nonlinear model function (nlsRespMod only)
##' @param pnames character vector of parameter names for the
##'        nonlinear model
##' @return an respModule object
mkRespMod <- function(fr, family = NULL, nlenv = NULL, nlmod = NULL)
{
    N <- n <- nrow(fr)
    if (!is.null(nlmod)) {
        nleta <- eval(nlmod, nlenv)
        grad <- attr(nleta, "gradient")
        if (is.null(grad))
            stop("At present a nonlinear model must return a gradient attribute")
        N <- n * ncol(grad)
    }
                                        # components of the model frame
    y <- model.response(fr)
    if(length(dim(y)) == 1) { # avoid problems with 1D arrays, but keep names
        nm <- rownames(y)
        dim(y) <- NULL
        if(!is.null(nm)) names(y) <- nm
    }
    weights <- model.weights(fr)
    if (is.null(weights)) weights <- rep.int(1, n)
    else if (any(weights < 0))
        stop(gettext("negative weights not allowed", domain = "R-Matrix"))
    offset <- model.offset(fr)
    if (is.null(offset)) offset <- numeric(N)
    if (length(offset) == 1) offset <- rep.int(offset, N)
    else if (length(offset) != N)
        stop(gettextf("number of offsets (%d) should be %d (s * n)",
                      length(offset), N), domain = "R-Matrix")
    ll <- list(weights = unname(weights), offset = unname(offset),
               wtres = numeric(n))
    if (!is.null(family)) {
        ll$y <- y                       # may get overwritten later
        rho <- new.env()
        rho$etastart <- model.extract(fr, "etastart")
        rho$mustart <- model.extract(fr, "mustart")
        rho$nobs <- n
        if (is.character(family))
            family <- get(family, mode = "function", envir = parent.frame(3))
        if (is.function(family)) family <- family()
        eval(family$initialize, rho)
        family$initialize <- NULL       # remove clutter from str output
        ll$mu <- unname(rho$mustart)
        lr <- as.list(rho)
        ll[names(lr)] <- lr             # may overwrite y, weights, etc.
        ll$weights <- unname(ll$weights)
        ll$y <- unname(ll$y)
        ll$eta <- family$linkfun(ll$mu)
        ll$sqrtrwt <- sqrt(ll$weights/family$variance(ll$mu))
        ll$sqrtXwt <- matrix(ll$sqrtrwt * family$mu.eta(ll$eta))
        ll$family <- family
        ll <- ll[intersect(names(ll), slotNames("glmRespMod"))]
        ll$n <- unname(rho$n)           # for the family$aic function
        ll$Class <- "glmRespMod"
    } else {
        ll$sqrtrwt <- sqrt(ll$weights)
        ll$y <- unname(as.numeric(y))
        ll$mu <- numeric(n)
        if (is.null(nlenv)) {
            ll$Class <- "respModule"
            ll$sqrtXwt <- matrix(ll$sqrtrwt)
        } else {
            ll$Class <- "nlsRespMod"
            ll$nlenv <- nlenv
	    ll$nlmod <- quote(nlmod)
            ll$sqrtXwt <- grad
            ll$pnames <- colnames(ll$sqrtXwt)
        }
    }
    do.call("new", ll)
}

glm4 <- function(formula, family, data, weights, subset,
                 na.action, start = NULL, etastart, mustart, offset,
		 sparse = FALSE, drop.unused.levels = FALSE, doFit = TRUE,
		 control = list(...),
                 ## all the following are currently ignored:
                 model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...) {
    call <- match.call()
    if (missing(family)) {
        family <- NULL
    } else {
        if(is.character(family))
            family <- get(family, mode = "function", envir = parent.frame())
        if(is.function(family)) family <- family()
        if(is.null(family$family)) {
            print(family)
            stop("'family' not recognized")
        }
    }
    ## extract x, y, etc from the model formula and frame
    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ans <- new("glpModel", call = call,
	       resp = mkRespMod(mf, family),
	       pred = as(model.Matrix(formula, mf, sparse = sparse,
				      drop.unused.levels=drop.unused.levels),
			 "predModule"))
    if (doFit)
	## TODO ? - make 'doFP' a function argument / control component:
	fitGlm4(ans, doFP = TRUE, control = control)
    else
	ans
}

fitGlm4 <- function(lp, doFP = TRUE, control = list()) {
### note that more than one iteration would need to update more than just 'coef'
    if(doFP && is(lp@resp, "glmRespMod"))
        lp@pred@coef <- glm.fp(lp)
    IRLS(lp, control)
}

##' A single step in the fixed-point algorithm for GLMs.
##'
##' In general we use an algorithm similar to the Gauss-Newton
##' algorithm for nonlinear least squares (except, of course, that it
##' allows for reweighting).  For some models, such as those using the
##' Gamma family with the inverse link the initial values of eta must
##' be non-zero.  This function calculates a single iteration of the
##' fixed-point algorithm used in stats::glm.fit to obtain suitable
##' starting estimates for the parameters.
##' @title Fixed-point iteration for a GLM
##' @param lp a linear predictor model.  The resp slot should inherit
##' from the glmRespMod class.
##' @return parameter vector
glm.fp <- function(lp) {
    stopifnot(is(lp, "glpModel"), is(rM <- lp@resp, "glmRespMod"))
    ff <- rM@family
    mu <- rM@mu
    vv <- ff$variance(mu)
    eta <- rM@eta
    muEta <- ff$mu.eta(eta)
    wts <- rM@weights
    z <- (eta - rM@offset) + (rM@y - rM@mu)/muEta
    good <- is.finite(vv) & vv > 0 & is.finite(z)
    stopifnot(any(good))
    w <- sqrt(wts * muEta * muEta /vv)[good]
    wM <- lp@pred@X[good,] * w
    as.vector(solve(crossprod(wM), crossprod(wM, z[good] * w)))
}

##'
##' @title
##' @param control  a (named) list {or vector; as.list(.)  must work}.
##' @param defaults a (named) list {or vector; as.list(.)  must work}.
##' @param rho typically an environment; in fact anything that "works" as third
##'    argument in  'assign(nm, val, rho)'>
##' @param nomatch.action string specifying what should happen when control()
##'    entries do not match any of the defaults.
##' @return none. Side effect: 'rho' will contain 'control' and 'defaults' entries.
##' @author Doug Bates (& Martin Maechler)
do.defaults <- function(control, defaults, rho,
			## by default stop() on mistyped control arguments:
			nomatch.action = c("stop", "warning", "none"))
{
    nomatch.action <- match.arg(nomatch.action)
					# Install the default values
    dnms <- names(defaults <- as.list(defaults))
    lapply(dnms, function(nm) assign(nm, defaults[[nm]], rho))
					# Match names of control arguments to defaults
    matched <- !is.na(mm <- pmatch(names(control <- as.list(control)), dnms))
    if(nomatch.action != "none" && any(!matched)) {
	msg <-
	    paste("The following control arguments did not match any default's names:",
		  paste(dQuote(names(control)[!matched]), collapse=", "), sep="\n   ")
	switch(nomatch.action,
	       "warning" = warning(msg, call.=FALSE, immediate.=TRUE),
	       "stop" = stop(msg, call.=FALSE))
    }
    if (any(matched)) {
	cc <- control[matched]
	names(cc) <- dnms[mm[matched]]
	lapply(names(cc),
	       function(nm) assign(nm, as(cc[[nm]], class(defaults[[nm]])), rho))
    }
    invisible()
}

IRLS <- function(mod, control) {
    stopifnot(is(mod, "glpModel"))
    respMod <- mod@resp
    predMod <- mod@pred
    ## localVariables("..."):
    MXITER <- warnOnly <- verbose <- quick <- TOL <- SMIN <- finalUpdate <- NULL
    do.defaults(control,
		list(MXITER = 200L, TOL = 0.0001, SMIN = 0.0001,
		     verbose = 0L,# integer: for verboseness levels
		     warnOnly = FALSE,
		     quick = TRUE, finalUpdate = FALSE),
                environment())
    cc <- predMod@coef
    respMod <- updateMu(respMod, as.vector(predMod@X %*% cc))
    iter <- nHalvings <- 0 ; DONE <- FALSE
    repeat {
	if((iter <- iter + 1) > MXITER) {
            msg <- paste("Number of iterations exceeded maximum MXITER =", MXITER)
            if(!warnOnly)
                stop(msg)
            ## else :
            warning(msg)
            cc <- cbase
            DONE <- TRUE
            break
        }
        cbase <- cc
        respMod <- updateWts(respMod)
        wrss0 <- sum(respMod@wtres^2)
        predMod <- reweightPred(predMod, respMod@sqrtXwt, respMod@wtres)
        incr <- solveCoef(predMod)
        convcrit <- sqrt(attr(incr, "sqrLen")/wrss0)
	if(verbose)
	    cat(sprintf("_%d_ convergence criterion: %5g\n",
			iter, convcrit))
        if(quick)## faster, but "loses" precision by not doing the "free" update:
            if (convcrit < TOL) break
        step <- 1
        repeat {
            cc <- as.vector(cbase + step * incr)
            respMod <- updateMu(respMod, as.vector(predMod@X %*% cc))
            wrss1 <- sum(respMod@wtres^2)
            if (verbose) {
		cat(sprintf("step = %.5f, new wrss = %.8g, Delta(wrss)= %g, coef =\n",
                            step, wrss1, wrss0 - wrss1))
                print(cc)
            }
	    if (wrss1 < wrss0) break
	    ## else
	    if ((step <- step/2) < SMIN) {
                msg <- "Minimum step factor 'SMIN' failed to reduce wrss"
		if(!warnOnly)
                    stop(msg)
                ## else :
                warning(msg)
                cc <- cbase
                DONE <- TRUE
                break
	    }
            ## no further step halving, if we are good enough anyway
	    if (DONE <- convcrit < TOL) break
            nHalvings <- nHalvings + 1
	}
        if(DONE || (!quick # check now
                    && convcrit < TOL))
            break
    }
    predMod@coef <- cc
    if(finalUpdate) {
	respMod <- updateWts(respMod)
	predMod <- reweightPred(predMod, respMod@sqrtXwt, respMod@wtres)
    }

    mod@ fitProps <- list(convcrit=convcrit, iter=iter, nHalvings=nHalvings)
    ## This is more portable than  new("glpModel", ....) as soon as
    ## the class contains extra slots (such as 'call'):
    mod@ resp <- respMod
    mod@ pred <- predMod
    mod
}

setMethod("formula", "Model", function(x, ...) x@call$formula)
setMethod("coef", "glpModel", function(object, ...)
      {
	  prd <- object@pred
	  structure(prd@coef,
		    names = colnames(prd@X))
      })
setMethod("fitted", "respModule", function(object, ...) object@mu)
setMethod("fitted", "glpModel", function(object, ...) {object <- object@resp; callGeneric(...)})

setMethod("residuals", "respModule",
          function(object, type = c("deviance", "pearson",
                           "working", "response", "partial"), ...)
      {
	  type <- match.arg(type)
          if (type %in% c("pearson", "deviance")) return(object@wtres)
          if (type %in% c("working", "response")) return(object@y - object@mu)
          stop(paste("residuals of type", sQuote(type), "not yet available"))
      })
setMethod("residuals", "glmRespMod",
          function(object, type = c("deviance", "pearson",
                           "working", "response", "partial"), ...)
      {
	  type <- match.arg(type)
          if (type == "pearson") return(object@wtres)

          fam <- object@family
	  mu <- object@mu
	  y <- object@y
          wts <- object@weights
          residuals <- y - mu
	  if (type == "response") return(residuals)
          if (type == "working") return(residuals/fam$mu.eta(object@eta))
          if (type == "deviance") {
              d.res <- sqrt(pmax(fam$dev.resids(y, mu, wts), 0))
              return(ifelse(y > mu, d.res, -d.res))
          }
          stop(paste("residuals of type", sQuote(type), "not yet available"))
      })
setMethod("residuals", "glpModel",
          function(object, type = c("deviance", "pearson",
                           "working", "response", "partial"), ...)
      {
          object <- object@resp
          callGeneric(...)
      })

setMethod("updateMu", signature(respM = "respModule", gamma = "numeric"),
	  function(respM, gamma, ...)
      {
	  respM@ wtres <- respM@sqrtrwt *
	      (respM@y - (respM@ mu <- respM@offset + gamma))
	  respM
      })

setMethod("updateMu", signature(respM = "glmRespMod", gamma = "numeric"),
          function(respM, gamma, ...)
      {
          respM@ mu <- respM@family$linkinv(respM@ eta <- respM@offset + gamma)
          respM@ wtres <- respM@sqrtrwt * (respM@y - respM@mu)
          respM
      })

setMethod("updateMu", signature(respM = "nlsRespMod", gamma = "numeric"),
          function(respM, gamma, ...)
      {
          ll <- as.data.frame(matrix(respM@offset + gamma,
                                     nrow = length(respM@y),
                                     dimnames = list(NULL, respM@pnames)))
          lapply(names(ll),
                 function(nm) assign(nm, ll[[nm]], envir = respM@nlenv))
          mm <- eval(respM@nlmod, respM@nlenv)
          respM@ wtres <- respM@sqrtrwt * (respM@y - (respM@ mu <- as.vector(mm)))
          respM@ sqrtXwt <- respM@sqrtrwt * attr(mm, "grad")
          respM
      })
setMethod("updateMu", signature(respM = "nglmRespMod", gamma = "numeric"),
	  function(respM, gamma, ...)
      {
	  .NotYetImplemented() ## FIXME
      })


## For models based on a Gaussian distribution (incl. "nlsRespMod")
## updateWts() has no effect:
setMethod("updateWts", signature(respM = "respModule"),
          function(respM, ...) respM)

setMethod("updateWts", signature(respM = "glmRespMod"),
          function(respM, ...)
      {
	  respM@ sqrtrwt   <- rtrwt <- sqrt(respM@weights/respM@family$variance(respM@mu))
	  respM@ sqrtXwt[] <- rtrwt * respM@family$mu.eta(respM@eta)
	  respM@ wtres	   <- rtrwt * (respM@y - respM@mu)
	  respM
      })

setMethod("reweightPred",
          signature(predM = "dPredModule", sqrtXwt = "matrix", wtres = "numeric"),
          function(predM, sqrtXwt, wtres, ...)
      {
          V <- as.vector(sqrtXwt) * predM@X
          s <- ncol(sqrtXwt)
          if (s > 1L)
              V <- Reduce("+", lapply(split(seq_len(nrow(V)), gl(s, nrow(sqrtXwt))),
                                      function(ind) V[ind,]))
          predM@Vtr <- as.vector(crossprod(V, wtres))
          predM@fac <- chol(crossprod(V))
          predM
      })

setMethod("reweightPred",
          signature(predM = "sPredModule", sqrtXwt = "matrix", wtres = "numeric"),
          function(predM, sqrtXwt, wtres, ...)
      {
          Vt <- crossprod(predM@X, Diagonal(x = as.vector(sqrtXwt)))
          s <- ncol(sqrtXwt)
          if (s > 1L)
              Vt <- Reduce("+", lapply(split(seq_len(ncol(Vt)), gl(s, nrow(sqrtXwt))),
                                      function(ind) Vt[, ind]))
          predM@Vtr <- as.vector(Vt %*% wtres)
          predM@fac <- update(predM@fac, Vt)
          predM
      })

setMethod("solveCoef", "dPredModule", function(predM, ...)
      {
          cc <- solve(t(predM@fac), predM@Vtr)
	  structure(as.vector(solve(predM@fac, cc)),
		    sqrLen = sum(as.vector(cc)^2))
      })

setMethod("solveCoef", "sPredModule", function(predM, ...)
      {
          ff <- predM@fac
          if (isLDL(ff)) stop("sparse factor must be LL, not LDL")
          cc <- solve(ff, solve(ff, predM@Vtr, system = "P"), system = "L")
	  structure(as.vector(solve(ff, solve(ff, cc, system = "Lt"), system = "Pt")),
		    sqrLen = sum(as.vector(cc)^2))
      })
