
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                 DESCRIPTION:
#  isPositiveDefinite        M  Checks if the matrix X is positive definite
#  makePositiveDefinite      M  Forces the matrix x to be positive definite
################################################################################


isPositiveDefinite <-
    function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Checks if the matrix x is positive definite

    # Arguments:
    #   x - a symmetric matrix or any other rectangular object
    #       describing a covariance matrix which can be converted into
    #       a matrix by the function 'as.matrix'.

    # FUNCTION:

    # Transform:
    x = as.matrix(x)

    # Check if matrix is positive definite:
    ans = .is.positive.definite(m = x)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.is.positive.definite <-
    function (m, tol, method = c("eigen", "chol"))
{
    # Author:
    #   Copyright 2003-05 Korbinian Strimmer
    #   Rank, condition, and positive definiteness of a matrix
    #   GNU General Public License, Version 2

    method = match.arg(method)
    if (!is.matrix(m)) {
        m = as.matrix(m)
    }
    if (method == "eigen") {
        eval = eigen(m, only.values = TRUE)$values
        if( missing(tol) ) {
            tol = max(dim(m))*max(abs(eval))*.Machine$double.eps
        }
        if (sum(eval > tol) == length(eval)) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    } else if (method == "chol") {
        val = try(chol(m), silent = TRUE)
        if (class(val) == "try-error") {
            return(FALSE)
        } else {
            return(TRUE)
        }
    }
}


# ------------------------------------------------------------------------------


makePositiveDefinite <-
    function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Forces the matrix x to be positive definite

    # Arguments:
    #   x - a symmetric matrix or any other rectangular object
    #       describing a covariance matrix which can be converted into
    #       a matrix by the function 'as.matrix'.

    # FUNCTION:

    # Make Positive Definite:
    ans = .make.positive.definite(m = x)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.make.positive.definite <-
    function(m, tol)
{
    # Author:
    #   Copyright 2003-05 Korbinian Strimmer
    #   Rank, condition, and positive definiteness of a matrix
    #   GNU General Public License, Version 2

    # Method by Higham 1988

    if (!is.matrix(m)) {
        m = as.matrix(m)
    }

    d = dim(m)[1]
    if ( dim(m)[2] != d ) {
        stop("Input matrix is not square!")
    }

    es = eigen(m)
    esv = es$values

    if (missing(tol)) {
        tol = d*max(abs(esv))*.Machine$double.eps
    }
    delta =  2*tol
        # factor two is just to make sure the resulting
        # matrix passes all numerical tests of positive definiteness

    tau = pmax(0, delta - esv)
    dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)

    # print(max(DA))
    # print(esv[1]/delta)

    return( m + dm )
}


################################################################################

