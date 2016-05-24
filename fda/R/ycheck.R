ycheck = function(y, n) {

#  check Y

if (is.vector(y)) y <- as.matrix(y)

if (!inherits(y, "matrix") && !inherits(y, "array"))
    stop("Y is not of class matrix or class array.")

ydim = dim(y)

if (ydim[1] != n) stop("Y is not the same length as ARGVALS.")

#  set number of curves and number of variables

ndim  = length(ydim)
if (ndim == 2) {
        ncurve = ydim[2]
        nvar   = 1
}
if (ndim == 3) {
        ncurve = ydim[2]
        nvar   = ydim[3]
}
if (ndim > 3) stop("Second argument must not have more than 3 dimensions")


return(list(y=y, ncurve=ncurve, nvar=nvar, ndim=ndim))

}


