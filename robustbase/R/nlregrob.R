#### nlrob.<meth>() functions for  high breakdown point  nlrob() methods

## concept (and original version) from lme4/R/lmer.R
getOptfun <- function(optimizer, needArgs = c("fn","par","lower","control"))
{
    if (((is.character(optimizer) && optimizer=="optimx") ||
         deparse(substitute(optimizer))=="optimx") &&
        !("package:optimx") %in% search())
        stop(shQuote("optimx")," package must be loaded in order to ",
             "use ",shQuote('optimizer="optimx"'))
    optfun <- if (is.character(optimizer))
	tryCatch(get(optimizer), error=function(e) NULL) else optimizer
    if (is.null(optfun))
        stop("couldn't find optimizer function ",optimizer )
    if (!is.function(optfun)) stop("non-function specified as optimizer")
    if (any(is.na(match(needArgs, names(formals(optfun))))))
	stop("optimizer function must use (at least) formal parameters ",
	     pasteK(sQuote(needArgs)))
    optfun
}

##' Utility for all nlrob.<meth>():  Find how and where to get parameter
##' names from, also check lower, upper, and replicate if needed.
##'
##' @param lower possibly unnamed numeric vector
##' @param upper as \code{lower}; both will be replicated to
##'    \code{length(pnames)} if that is specified and longer.
##' @param pnames DEPRECATED possibly missing character vector
##' @param var.nms character vector of which 'pnames' must be a subset of.
##' @param envir \code{\link{environment}: the function possibly assigns
##'    "lower", "upper" in the environment \code{envir}.
.fixupArgs <- function(lower, upper, pnames, var.nms, envir) {
    if(missing(pnames)) {
        if(is.null(pnames <- names(lower))) pnames <- names(upper)
        if(is.null(pnames))
            stop("Either specify 'pnames' or provide 'upper' or 'lower' with names()")
    } else if (!is.character(pnames))
        stop("'pnames' must be a character vector")
    else
        warning("specifying 'pnames' is deprecated; rather 'lower' or 'upper' should have names()")
    if(any(is.na(match(pnames, var.nms))))
        stop("parameter names must appear in 'formula'")
    npar <- length(pnames)
    if (npar > 1 && length(lower) == 1)
        envir$lower <- rep.int(lower, npar)
    else if (length(lower) != npar)
        stop(gettextf("lower must be either of length %d, or length 1", npar))
    if (npar > 1 && length(upper) == 1)
        envir$upper <- rep.int(upper, npar)
    else if (length(upper) != npar)
        stop(gettextf("upper must be either of length %d, or length 1", npar))
    stopifnot(is.numeric(lower), is.numeric(upper), lower <= upper)
    pnames
}

nlrob.MM <-
    function(formula, data, pnames, lower, upper, tol = 1e-6,
	     psi = c("bisquare", "lqq", "optimal", "hampel"),
	     init = c("S", "lts"),
             ctrl = nlrob.control("MM", psi=psi, init=init, fnscale=NULL,
                 tuning.chi.scale = .psi.conv.cc(psi, .Mchi.tuning.defaults[[psi]]),
                 tuning.psi.M     = .psi.conv.cc(psi, .Mpsi.tuning.defaults[[psi]]),
                 optim.control = list(), optArgs = list(...)),
             ...)
{
    ctrl.exp <- substitute(ctrl)
    if(missing(ctrl)) {
        init <- match.arg(init)
        psi  <- match.arg(psi)
        force(ctrl) #
    } else {
        init <- ctrl$ init
        psi  <- ctrl$ psi
    }
    c1 <- ctrl$tuning.chi.scale
    c2 <- ctrl$tuning.psi.M
    if(is.character(ctrl$optimizer)) {
### TODO
    } else if(is.function(ctrl$optimizer)) {
### TODO
    } else
	stop(gettextf("'%s' must be character string or function, but is \"%s\"",
		      "ctrl$optimizer", class(ctrl$optimizer)), domain=NA)


    ## Preliminary psi-specific checks / computations:
    switch(psi,
           "lqq" = { # lqqMax = rho(Inf), used in rho.inv() *and* 'constant':
               c12 <- c1[1]+c1[2]
               lqqMax <- (c1[1]*c1[3] - 2*c12)/(1-c1[3]) + c12})

    rho1 <- function(t) Mchi(t, c1, psi)
    rho2 <- function(t) Mchi(t, c2, psi)

    rho.inv <- switch(psi, "bisquare" = function(y) {
        ## Find  x := u^2  which solves cubic eq. 3*x - 3*x^2 + x^3 = y
        ##      <==> (x-1)^3 + 1 = y  <==> (1-x)^3 = 1-y  <==> x = 1 - (1-y)^(1/3)
        ##      (where we assume 0 <= y <= 1, i.e, y-1 < 0)
        c1 * sqrt(1 - (1 - y)^(1/3))
    }, "lqq" = function(y) {
        uniroot( function(x) rho1(x) - y, lower = 0, upper = lqqMax )$root
    }, "optimal" = function(y) {
        ## Salibian-Barrera, Matias, Willems, Gert, and Zamar, Ruben (2008).
        ## The fast-tau estimator for regression.
        ## Journal of Computational and Graphical Statistics 17, 659-682.
        sqrt(y/1.38) * c1 * 3
    }, "hampel" = function(y) {
        C <- MrhoInf(c1, psi)
        a <- c1[1]; b <- c1[2]; r <- c1[3]
        if (y <= a/C)
            sqrt(2*C*y)
        else if (y <= (2*b - a)/C)
            0.5*a + C/a*y
        else r + sqrt( r^2 - ( (r - b)*(2*C/a*y + (b - a)) - b*r ) )
    }, stop(gettextf("Psi function '%s' not supported yet", psi)))

    M_scale <- function(sigma, u) sum( rho1(u/sigma) )/nobs - 0.5

    objective.initial <-
        switch(init,
               "lts" = function(par) { ## and (h, formula, data, pnames)
                   y.hat <- eval( formula[[3L]], c(data, setNames(par, pnames)) )
                   sum(sort.int( (y - y.hat)^2, partial = h )[1:h])
               },
               "S" = function(par) { ## and (constant, formula, data, pnames)
                   y.hat <- eval( formula[[3L]], c(data, setNames(par, pnames)) )
                   res <- y - y.hat
                   ## Rousseeuw, Peter J., and Leroy, Annick M. (1987).
                   ## Robust Regression & Outlier Detection.
                   ## John Wiley & Sons, New York, p. 137.
                   med_abs_res <- median(abs(res))
                   uniroot(M_scale,
                           lower = constant[1L] * med_abs_res,
                           upper = constant[2L] * med_abs_res,
                           u = res )$ root ## == 'sigma'
	       }, stop(gettextf("Initialization 'init = \"%s\"' not supported (yet)",
				init)))

    objective.M <- function(par, sigma) {
        y.hat <- eval( formula[[3L]], c(data, setNames(par, pnames)) )
        sum(rho2( (y - y.hat)/sigma ))
    }
    ## => psi(.) / wgt(.) for robustness weights is
    ##  Mpsi(x, c2, psi) or Mwgt(x, c2, psi)

    formula <- as.formula(formula)
    dataName <- substitute(data)
    varNames <- all.vars(formula)
    obsNames <- rownames(data <- as.data.frame(data))
    data <- as.list(data)# to be used as such
    if (length(formula) == 2L) { ## as nls
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 0
    }

    npar <- length(pnames <- .fixupArgs(lower, upper, pnames, varNames, environment()))
    ##			      ^^^^^^^^^ -> possibly changes (lower, upper) in envir.
    y <- eval(formula[[2L]], data)
    nobs <- length(y)
    stopifnot(nobs >= npar)
    if (is.null(fnscale <- ctrl$ fnscale))
        fnscale <- sum((y - mean(y))^2)
    ctrl$fnscale <- NULL # remove it there
    stopifnot(is.numeric(fnscale), fnscale > 0)
    ## is used in M_scale() in any case, and in init-estim. if "S"
    constant <- c(
        switch(psi, bisquare = 1/c1,
               lqq = 1/lqqMax,
               optimal = 1/c1 * 1/3,
               hampel = 1/c1[3]),
        if(nobs %% 2) 2/rho.inv(2/(nobs+2)) else 1/rho.inv(1/(nobs+1)))
    switch(init, lts = h <- (nobs + npar + 1)%/%2)

    ## FIXME: "optimizer":
    initial <- do.call(JDEoptim, c(list(lower, upper, objective.initial,
					tol=tol, fnscale=fnscale), ctrl$optArgs))
    names(initial$par) <- pnames
    res <- y - eval( formula[[3L]], c(data, initial$par) )

    med_abs_res <- median(abs(res))
    sigma <- uniroot( M_scale,
                      lower = constant[1L] * med_abs_res,
                      upper = constant[2L] * med_abs_res,
                      u = res )$root

    M <- optim(initial$par, objective.M, sigma = sigma,
               method = "L-BFGS-B", lower = lower, upper = upper,
               control = c(list(fnscale = initial$value, parscale = initial$par),
			   ctrl$optim.control), hessian = TRUE)
    ## 'hessian': experimental - FIXME: eliminate if unused
    coef <- setNames(M$par, pnames)
    status <-
	if (M$convergence == 0) "converged"
	else if (M$convergence == 1)
	    "maximum number of iterations reached without convergence"
	else M$message
    fit <- eval( formula[[3L]], c(data, coef) )
    names(fit) <- obsNames
    structure(list(call = match.call(), formula=formula, nobs=nobs,
                   coefficients = coef,
                   fitted.values = fit,
                   residuals = y - fit,
                   crit = M$value,
                   initial = initial,
                   Scale = sigma,
                   status = status, counts = M$counts, data = dataName,
                   hessian = M$hessian, ctrl=ctrl),
              class = "nlrob")
} ## nlrob.MM


nlrob.tau <- function(formula, data, pnames, lower, upper, tol = 1e-6,
		      psi = c("bisquare", "optimal"),
		      ctrl = nlrob.control("tau", psi=psi, fnscale=NULL,
			  tuning.chi.scale = NULL, tuning.chi.tau = NULL,
			  optArgs = list(...)),
		      ...)
{
    if(missing(ctrl)) {
	psi <- match.arg(psi)
	force(ctrl) #
    } else {
	psi  <- ctrl$ psi
    }
    if(is.null(.chi.s <- ctrl$tuning.chi.scale))
	.chi.s <- switch(psi, bisquare = list(b = 0.20, cc = 1.55),
			 optimal = list(b = 0.5, cc = 0.405))
    if(is.null(.chi.t <- ctrl$tuning.chi.tau))
	.chi.t <- switch(psi, bisquare = list(b = 0.46, cc = 6.04),
			 optimal = list(b = 0.128, cc = 1.060))

    b1 <- .chi.s$b
    c1 <- .chi.s$cc
    b2 <- .chi.t$b
    c2 <- .chi.t$cc

    ## Preliminary psi-specific checks / computations:
    switch(psi, "bisquare" = {
	b1 <- b1/MrhoInf(c1, psi)
	b2 <- b2/MrhoInf(c2, psi)
    })

    rho1 <- function(t) Mchi(t, c1, psi)
    rho2 <- function(t) Mchi(t, c2, psi)

    rho.inv <- switch(psi, "bisquare" = function(y) {
	c1 * sqrt(1 - (1 - y)^(1/3))
    }, "optimal" = function(y) {
        ## Salibian-Barrera, Matias, Willems, Gert, and Zamar, Ruben (2008).
        ## The fast-tau estimator for regression.
        ## Journal of Computational and Graphical Statistics 17, 659-682.
        sqrt(y/1.38) * c1 * 3
    })

    M_scale <- function(sigma, u) sum( rho1(u/sigma) )/nobs - b1
    tau_scale2 <- function(u, sigma) sigma^2 * 1/b2*sum( rho2(u/sigma) )/nobs

    objective <- function(par) {
        fit <- eval( formula[[3L]], c(data, setNames(par, pnames)) )
        res <- y - fit
        ## Rousseeuw, Peter J., and Leroy, Annick M. (1987).
        ## Robust Regression & Outlier Detection.
        ## John Wiley & Sons, New York, p. 137.
        med_abs_res <- median(abs(res))
        sigma <- uniroot( M_scale,
                          lower = constant[1L] * med_abs_res,
                          upper = constant[2L] * med_abs_res,
                          u = res )$root
        tau_scale2(res, sigma)
    }

    formula <- as.formula(formula)
    dataName <- substitute(data)
    varNames <- all.vars(formula)
    obsNames <- rownames(data <- as.data.frame(data))
    data <- as.list(data)# to be used as such
    if (length(formula) == 2L) { ## as nls
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 0
    }

    npar <- length(pnames <- .fixupArgs(lower, upper, pnames, varNames, environment()))
    ##			      ^^^^^^^^^ -> possibly changes (lower, upper) in envir.
    y <- eval(formula[[2L]], data)
    nobs <- length(y)
    stopifnot(nobs >= npar)
    if (is.null(fnscale <- ctrl$ fnscale))
        fnscale <- mean((y - mean(y))^2)
    ctrl$fnscale <- NULL # remove it there
    stopifnot(is.numeric(fnscale), fnscale > 0)
    constant <- c(
        switch(psi,
               bisquare = 1/c1,
               optimal = 1/c1 * 1/3),
        if (nobs %% 2) 2/rho.inv(2/(nobs+2)) else 1/rho.inv(1/(nobs+1)))

    optRes <- do.call(JDEoptim, c(list(lower, upper, objective, tol=tol, fnscale=fnscale),
				  ctrl$optArgs))
    iter <- optRes$iter
    status <-
        if (optRes$convergence == 0)
            "converged"
        else paste("failed to converge in", iter, "steps")
    coef <- setNames(optRes$par, pnames)
    fit <- eval( formula[[3L]], c(data, coef) )
    names(fit) <- obsNames
    structure(list(call = match.call(), formula=formula, nobs=nobs,
                   coefficients = coef,
                   fitted.values = fit,
                   residuals = y - fit,
                   crit = optRes$value,
                   Scale = sqrt(optRes$value),
                   status = status, iter = iter, data = dataName, ctrl=ctrl),
              class = "nlrob")
} ## nlrob.tau


nlrob.CM <- function(formula, data, pnames, lower, upper, tol = 1e-6,
		     psi = c("bisquare", "lqq", "welsh", "optimal", "hampel", "ggw"),
                     ctrl = nlrob.control("CM", psi=psi, fnscale=NULL,
                         tuning.chi = NULL, optArgs = list(...)),
                     ...)
{
    if(missing(ctrl)) {
        psi <- match.arg(psi)
        force(ctrl) #
    } else {
        psi  <- ctrl$ psi
    }
    if (is.null(t.chi <- ctrl$tuning.chi))
	t.chi <- switch(psi, bisquare = list(b = 0.5, cc = 1, c = 4.835),
			stop("unable to find constants for psi function"))
    ## FIXME:
    b  <- t.chi$b  ## b = epsilon (in paper) = fraction of outlier ~= breakdown
    cc <- t.chi$cc ## cc = k; make
    c  <- t.chi$c  ## c = the factor in objective   c*rho(.) - log(sigma)

    rho <- function(t) Mchi(t, cc, psi)
    M_scale <- function(sigma, u) sum( rho(u/sigma) )/nobs - b

    formula <- as.formula(formula)
    dataName <- substitute(data)
    varNames <- all.vars(formula)
    obsNames <- rownames(data <- as.data.frame(data))
    data <- as.list(data)# to be used as such
    if (length(formula) == 2L) { ## as nls
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 0
    }

    npar <- length(pnames <- .fixupArgs(lower,upper,pnames, c(varNames,"sigma"),environment()))
    ##			      ^^^^^^^^^ -> possibly changes (lower, upper) in envir.
    if ("sigma" %in% pnames) {
	if ("sigma" %in% varNames || "sigma" %in% names(data))
	    stop("As \"sigma\" is in 'pnames', do not use it as variable or parameter name in 'formula'")
	stopifnot(lower[pnames == "sigma"] >= 0)

	objective <- function(par) {
	    par <- setNames(par, pnames)
	    fit <- eval( formula[[3L]], c(data, par) )
	    sigma <- par[["sigma"]]
	    c * sum(rho( (y - fit)/sigma ))/nobs + log(sigma)
	}
	con <- function(par) {
	    par <- setNames(par, pnames)
	    fit <- eval( formula[[3L]], c(data, par) )
	    M_scale(par[["sigma"]], y - fit)
	}
    } else { ## hmm, this case *really* is not CM properly
	objective <- function(par) {
	    fit <- eval( formula[[3L]], c(data, setNames(par, pnames)) )
	    resid <- y - fit
	    sigma <- mad(resid)
	    c * sum(rho( resid/sigma ))/nobs + log(sigma)
	}
	con <- NULL
    }

    y <- eval(formula[[2L]], data)
    nobs <- length(y)
    stopifnot(nobs >= npar)
    if (is.null(fnscale <- ctrl$ fnscale))
        fnscale <- mean((y - mean(y))^2)
    ctrl$fnscale <- NULL # remove it there
    stopifnot(is.numeric(fnscale), fnscale > 0)

    optRes <- do.call(JDEoptim, c(list(lower, upper, objective, constr=con,
				       tol=tol, fnscale=fnscale),
				  ctrl$optArgs))
    iter <- optRes$iter
    status <- if (optRes$convergence == 0)
        "converged"
    else paste("failed to converge in", iter, "steps")
    coef <- setNames(optRes$par, pnames)
    fit <- eval( formula[[3L]], c(data, coef) )
    names(fit) <- obsNames
    structure(list(call = match.call(), formula=formula, nobs=nobs,
                   coefficients = coef,
                   fitted.values = fit,
                   residuals = y - fit,
                   crit = optRes$value,
                   status = status, iter = iter, data = dataName, ctrl=ctrl),
              class = "nlrob")
} ## nlrob.CM


nlrob.mtl <- function(formula, data, pnames, lower, upper, tol = 1e-6,
                      ctrl = nlrob.control("mtl", cutoff = 2.5, optArgs = list(...)),
                      ...)
{
    stopifnot(is.numeric(cutoff <- ctrl[["cutoff"]]), length(cutoff) >= 1)
    trim <- function(t) { # t = residuals Res, or  Res / sigma
        t <- sort.int(t)
        i <- which(t >= cutoff)
        h <- if (length(i) > 0)
            max(hlow, floor(min( (i - 1)/(2*pnorm(t[i]) - 1) ))) else nobs
        list(h = h, t = t)
    }

    formula <- as.formula(formula)
    dataName <- substitute(data)
    varNames <- all.vars(formula)
    obsNames <- rownames(data <- as.data.frame(data))
    data <- as.list(data)# to be used as such
    if (length(formula) == 2L) { ## as nls
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 0
    }

    npar <- length(pnames <- .fixupArgs(lower,upper,pnames, c(varNames,"sigma"),environment()))
    ##			      ^^^^^^^^^ -> possibly changes (lower, upper) in envir.
    constant <- log(2*pi)
    if ("sigma" %in% pnames) {
	if ("sigma" %in% varNames || "sigma" %in% names(data))
	    stop("As \"sigma\" is in 'pnames', do not use it as variable or parameter name in 'formula'")
	stopifnot(lower[pnames == "sigma"] >= 0)
	objective <- function(par) {
	    par <- setNames(par, pnames)
	    fit <- eval( formula[[3L]], c(data, par) )
	    sigma <- par[["sigma"]]
	    tp <- trim( abs( (y - fit)/sigma ) )
	    h <- tp$h
	    h*(constant + 2*log(sigma)) + sum(tp$t[1L:h]^2)
	}
    } else { ## hmm... this is not really MTL
	objective <- function(par) {
	    fit <- eval( formula[[3L]], c(data, setNames(par, pnames)) )
	    resid <- y - fit
	    sigma <- mad(resid)
	    tp <- trim( abs(resid/sigma) )
	    h <- tp$h
	    h*(constant + 2*log(sigma)) + sum(tp$t[1L:h]^2)
	}
    }

    y <- eval(formula[[2L]], data)
    nobs <- length(y)
    stopifnot(nobs >= npar)
    if (is.null(fnscale <- ctrl$ fnscale))
        fnscale <- sum((y - mean(y))^2)
    ctrl$fnscale <- NULL # remove it there
    stopifnot(is.numeric(fnscale), fnscale > 0)
    hlow <- (nobs + npar + 1)%/%2

    optRes <- do.call(JDEoptim, c(list(lower, upper, objective, tol=tol, fnscale=fnscale),
				  ctrl$optArgs))
    coef <- setNames(optRes$par, pnames)
    crit <- optRes$value
    iter <- optRes$iter
    status <- if (optRes$convergence == 0)
        "converged"
    else paste("failed to converge in", iter, "steps")
    fit <- eval( formula[[3L]], c(data, coef) )
    names(fit) <- obsNames
    resid <- y - fit
    quan <-
        trim( resid/(if ("sigma" %in% pnames) coef["sigma"] else mad(resid)) )$h

    structure(list(call = match.call(), formula=formula, nobs=nobs,
                   coefficients = coef,
                   fitted.values = fit,
                   residuals = resid,
                   crit = crit,
                   quan = quan,
                   status = status, iter = iter, data = dataName, ctrl = ctrl),
              class = "nlrob")
} ## nlrob.mtl

nlrob.control <- function(method,
                          psi = c("bisquare", "lqq", "welsh", "optimal", "hampel", "ggw"),
                          init = c("S", "lts"),
                          optimizer = "JDEoptim", optArgs  = list(),
                          ...)
{
  psi <- match.arg(psi)
  init <- match.arg(init)
  dots <- list(...)
  argNms <- names(dots)
  ##' argument or default -> return list of length 1
  a. <- function(nm,def) {
    L <- list( if(nm %in% argNms) dots[[nm]] else def )
    names(L) <- nm
    L
  }
  switch(method,
         "M" = {
             list(method = method) # not yet used
         },
         "MM" = {
             c(list(method = method, init = init, psi = psi),
               a.("fnscale", NULL),
               a.("tuning.chi.scale", .psi.conv.cc(psi, .Mchi.tuning.defaults[[psi]])),
               a.("tuning.psi.M",     .psi.conv.cc(psi, .Mpsi.tuning.defaults[[psi]])),
               a.("optim.control", list()),
               list(optimizer = optimizer, optArgs = optArgs))
         },
         "tau" = {
             c(list(method = method, psi = psi),
               a.("fnscale", NULL),
               a.("tuning.chi.scale", NULL),
               a.("tuning.chi.tau", NULL),
               list(optimizer = optimizer, optArgs = optArgs))
         },
         "CM" = {
             c(list(method = method, psi = psi),
               a.("fnscale", NULL),
               a.("tuning.chi", NULL),
               list(optimizer = optimizer, optArgs = optArgs))
         },
         "mtl" = {
             c(list(method = method),
               a.("fnscale", NULL),
               a.("cutoff", 2.5),
               list(optimizer = optimizer, optArgs = optArgs))
         },
         stop("Method ", method, "not correctly supported yet"))
}
