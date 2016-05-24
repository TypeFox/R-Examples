sqsvmpath <- function(x, y, nlam, flmin, ulam, isd, 
    eps, dfmax, pmax, jd, pf, pf2, maxit, lam2, nobs, nvars, vnames) {
    #################################################################################
    #data setup
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    if (!all(y %in% c(-1, 1))) 
        stop("y should be a factor with two levels")
    #################################################################################
    # call Fortran core
    fit <- .Fortran("sqsvmlassoNET", lam2, nobs, nvars, as.double(x), 
        as.double(y), jd, pf, pf2, dfmax, pmax, nlam, flmin, ulam, 
        eps, isd, maxit, nalam = integer(1), b0 = double(nlam), 
        beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
        alam = double(nlam), npass = integer(1), jerr = integer(1), 
        PACKAGE = "gcdnet")
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) <- c("sqsvmpath")
    outlist
} 
