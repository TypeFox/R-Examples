##
##  f m i n s e a r c h . R
##


fminsearch <- function(f, x0, ..., minimize = TRUE, dfree = TRUE,
                       maxiter = 1000, tol = .Machine$double.eps^(2/3)) {
    if (!is.numeric(x0))
        stop("Argument 'x0' must be a numeric vector.")

    scl <- if(minimize) 1 else -1
    fun <- match.fun(f)
    f <- function(x) scl * fun(x, ...)

    if (dfree) {
        opt <- nelder_mead(x0, f, tol = tol, maxfeval = 10*maxiter)
    } else {
        opt <- fletcher_powell(x0, f, maxiter = maxiter, tol = tol)
    }

    xopt <- opt$xmin; fopt <- opt$fmin
    if (! minimize) fopt <- -fopt
    return(list(xval = xopt, fval = fopt, niter = opt$niter))
}
