# Poisson Likelihood

sSolver <- function(z, a, extra) {
    x <- z
    if (is.null(extra$y)) stop("sSolver needs the additional argument y!")
    z <- extra$y
    fobj <- function(x) sum(x-z*log(x))       #target value
    gobj <- function(x) 1-z/x                 #gradient
    options(warn = -1)
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}
