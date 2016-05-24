#### Check tools -- notably for Rmpfr
#### ================================
## to be used as
##	source(system.file("check-tools.R", package="Rmpfr"), keep.source=FALSE)

## Get all the (non-Matrix specific) Matrix package check tools:
##	(needs  'Suggests: Matrix' in ../DESCRIPTION )
source(system.file("test-tools-1.R", package="Matrix"), keep.source=FALSE)
##    MM = ~/R/Pkgs/Matrix/inst/test-tools-1.R

### ------- Part I --  do not need 'Rmpfr'

`%=N=%` <- function(x,y) (x == y) | (is.na(x) & is.na(y))

all.eq.finite <- function(x,y, ...) {
    ## x = 'target'   y = 'current'
    if(any(is.finite(y[!(fx <- is.finite(x))])))
	return("current has finite values where target has not")
    if(any(is.finite(x[!(fy <- is.finite(y))])))
	return("target has finite values where current has not")
    ## now they have finite values at the same locations
    all.equal(x[fx], y[fy], ...)
}

### ------- Part II --  do not make sense or work outside of 'Rmpfr' :

all.EQ <- function(x,y, tolerance = 2^-98, ...) # very small tol. for MPFR
    all.equal.finite(x, y, tolerance=tolerance, ...)

