alspath <- function(x, y, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, 
                    pf, pf2, maxit, lam2, tau, nobs, nvars, vnames) {
    #################################################################################
    # data setup
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    if (tau <= 0 || tau >= 1) stop("tau must be in (0,1)")
	  tau <- as.double(tau)
    #################################################################################
    # call Fortran core
    fit <- .Fortran("alslassoNET", tau, lam2, nobs, nvars, x, 
        y, jd, pf, pf2, dfmax, pmax, nlam, flmin, ulam, 
        eps, isd, intr, maxit, nalam = integer(1), b0 = double(nlam), 
        beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
        alam = double(nlam), npass = integer(1), jerr = integer(1), 
        PACKAGE = "SALES")
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) <- c("alspath")
    outlist
} 
