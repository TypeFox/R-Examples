#### varNPreg.R : Nonparametric Variance Estimator
####
#### S/R interface to the resest() Fortran subroutine

#### Copyright © Martin Maechler (2001).
#### This software is distributed under the terms of the GNU GENERAL
#### PUBLIC LICENSE Version 2, June 1991, see the COPYING file from R,
#### or http://www.gnu.org/copyleft/gpl.html

varNPreg <- function(x,y)
{
    ## Purpose: Nonparametric Leave-1-out Residuals and Variance Estimator
    ##	in the model   y[i] = mu(x(i)) + E[i] ,  E[i] ~ (0, sigma^2), i.i.d

    ## Author: Martin Maechler, Date:  9 Jul 2001, 14:47
    if(2 >= (n <- length(x))) stop("n := length(x)  must be at least 3")
    if(is.unsorted(x)) stop("'x' must be ordered increasingly")
    if(n != length(y)) stop("'x' and 'y' must have same length")
    .Fortran(resest,
             as.double(x), as.double(y), n,
             res = double(n),
             snr = double(1),
	     sigma2 = double(1))[4:6]
}

varest <- function(x,y) {
    warning("Use  'varNPreg' instead of 'varest';\n",
	    "The latter has been an accidental identical copy of the former")
    varNPreg(x,y)
}
