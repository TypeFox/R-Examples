cpalspath <- function(x, y, w, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, 
                      pfmean, pf2mean, pfscale, pf2scale, maxit, lam2,
                      tau, nobs, nvars, vnames) {
    #################################################################################
    #data setup
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    if (tau <= 0 || tau >= 1 || tau == 0.5) 
        stop("tau must be in (0,1)\\{0.5}")
	tau <- as.double(tau)
    if (w <= 0) stop("Weight must be positive.")
    w <- as.double(w)
    #################################################################################
    # call Fortran core
    fit <- .Fortran("cpalslassoNET", w, tau, lam2, nobs, nvars, x, y, jd, pfmean, 
        pfscale, pf2mean, pf2scale, dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, 
        maxit, nalam = integer(1), b0 = double(nlam), beta = double(pmax * nlam), 
        ibeta = integer(pmax), nbeta = integer(nlam), t0 = double(nlam), 
        theta = double(pmax * nlam), itheta = integer(pmax), ntheta = integer(nlam), 
        alam = double(nlam), npass = integer(1), jerr = integer(1), PACKAGE = "SALES")
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) <- c("cpalspath")
    outlist
} 
