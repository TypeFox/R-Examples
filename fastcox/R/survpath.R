######################################################################
## These functions are minor modifications or directly copied from the 
## glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). 
##        Regularization Paths for Generalized Linear Models via 
##        Coordinate Descent. 
##        Journal of Statistical Software, 33(1), 1-22. 
##        URL http://www.jstatsoft.org/v33/i01/.
## The reason they are copied here is because they are internal functions
## and hence are not exported into the global environment.
## The original comments and header are preserved.

survpath <- function(x, y, d, nlam, flmin, ulam, isd, eps, dfmax, 
    pmax, jd, pf, maxit, alpha, nobs, nvars, vnames) {
    #################################################################################
    #data
    #   setup
    o <- order(d, decreasing = T)
    oo <- o[order(y[o])]
    x <- x[oo, ]
    y <- y[oo]
    d <- d[oo]
    rs <- which(d == 1)
    nrs <- length(rs)
    rs <- as.integer(rs)
    nrs <- as.integer(nrs)
    #################################################################################
    # call
    #   Fortran
    #   core
    fit <- .Fortran("coxlassoNET", rs, nrs, alpha, nobs, nvars, as.double(x), 
        jd, pf, dfmax, pmax, nlam, flmin, ulam, eps, isd, maxit, nalam = integer(1), 
        beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
        alam = double(nlam), npass = integer(1), jerr = integer(1), PACKAGE = "fastcox")
    #################################################################################
    #
    #   output
    errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
    switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = warning(errmsg$msg, 
        call. = FALSE))
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    class(outlist) <- c("survpath")
    outlist
}


 
