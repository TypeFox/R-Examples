tweediegrppath <- function(bn, bs, ix, iy, gamma, x, y, rho, vt, alpha, nlam, flmin, ulam, isd, eps, dfmax, pmax, pf, maxit, nobs, nvars, vnames, group) {
    #################################################################################
    # call Fortran core
    gamma <- 0.25 * gamma/nobs # not used
    gamma <- as.double(gamma) # not used
    fit <- .Fortran("tweediegrpnet", as.double(alpha), as.double(rho), as.double(vt), bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
        as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, as.integer(isd), maxit, nalam = integer(1), 
        b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
        nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
    outlist
}

