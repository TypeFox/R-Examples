olsmod <-
function (X, y, ind, tind, n, k, t, nT, w, w2, coef0 = rep(0, dim(X)[[2]]),
    hess = FALSE, trace = trace, x.tol = 1.5e-18, rel.tol = 1e-15,
    ...)
{

    ## extensive function rewriting, Giovanni Millo 29/09/2010
    ## structure:
    ## a) specific part
    ## - set names, bounds and initial values for parms
    ## - define building blocks for likelihood and GLS as functions of parms
    ## - define likelihood
    ## b) generic part(independent from ll.c() and #parms)
    ## - fetch covariance parms from max lik
    ## - calc last GLS step
    ## - fetch betas
    ## - calc final covariances
    ## - make list of results

    ## this function just for compatibility:
    ## ML estimation of OLS model; produces (hopefully same)
    ## results as lm() and logLik.lm but as a 'splm' object

    ## good for any GLS estimation by ML, just supply the
    ## right sigma.1() function

    ## set names for final parms vectors
    nam.beta <- dimnames(X)[[2]]
    #nam.errcomp <- c("psi")

    ## initialize values for optimizer
    ## (NULL is currently passed as coef0 if ols and lag=F)
    coef0 = rep(0, dim(X)[[2]])
    myparms0 <- coef0
    ## set bounds for optimizer
    lower.bounds <- -Inf
    upper.bounds <- Inf

    ## modules for likelihood
    ## (none)

    ## likelihood function, both steps included
    ll.c <- function(betas, y, X, n, t, w, w2) {
        ## get e, s2e as function of betas
        e <- y - X %*% betas                # lag-specific line (Ay for y)
        s2e <- crossprod(e)/(n*t)
        ## calc ll
        tre <- -n * t/2 * log(s2e)
        quattro <- -1/(2 * s2e) * crossprod(e)
        const <- -(n * t)/2 * log(2 * pi)
        ll.c <- const + tre + quattro
        ## invert sign for minimization
        llc <- -ll.c
    }


    ## GLS step function suppressed

    ## max likelihood
    optimum <- nlminb(start = myparms0, objective = ll.c,
                      gradient = NULL, hessian = NULL,
                      y = y, X = X, n = n, t = t, w = w, w2 = w2,
                      scale = 1, control = list(x.tol = x.tol,
                                 rel.tol = rel.tol, trace = trace),
                      lower = lower.bounds, upper = upper.bounds)

    ## log likelihood at optimum (notice inverted sign)
    myll <- -optimum$objective
    ## retrieve optimal parms
    betas <- optimum$par

    ## one last GLS step at optimal vcov parms suppressed

    ## final vcov(beta)
    e <- y - X %*% betas                # lag-specific line (Ay for y)
    s2e <- crossprod(e)/(n*t)
    covB <- as.numeric(s2e) * solve(crossprod(X))

    ## final vcov(errcomp)
    covAR <- NULL                                  # ols.errors-specific
    covPRL <- NULL                                 # ols.errors-specific

    ## final parms
    #betas ok
    arcoef <- NULL                                 # ols.errors-specific line
    errcomp <- NULL                                # ols.errors-specific
    names(betas) <- nam.beta
    #names(arcoef) <- "psi"                        # lag-specific line
    #names(errcomp) <- nam.errcomp[which(nam.errcomp!="psi")]

    dimnames(covB) <- list(nam.beta, nam.beta)
    #dimnames(covAR) <- list(names(arcoef), names(arcoef))
    #dimnames(covPRL) <- list(names(errcomp), names(errcomp))

    ## result
    RES <- list(betas = betas, arcoef=arcoef, errcomp = errcomp,
                covB = covB, covAR=covAR, covPRL = covPRL, ll = myll,
                sigma2 = s2e)

    return(RES)
}
