
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


################################################################################
# FUNCTION:                 DESCRIPTION:
#  norm2                      Returns the norm2 of a matrix
################################################################################


# IMPORTANT NOTICE:
#   The function fBasics::norm() has become obsolete, use instead: 
#   base::norm(x, type) - Note, the arguments are different.
#   The original Rmetrics function is still available as fBasics:norm2()


norm2 <-
function(x, p = 2)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the spectral norm of a matrix

    # Details:
    #   http://mathworld.wolfram.com/MatrixNorm.html:
    #   For p = 1
    #       The maximum absolute column sum norm |A|_1 is defined
    #       as the maximum of the sum of the absolute valued elements
    #       of columns of the matrix.
    #   For p = 2:
    #       The spectral |A|_2 norm is "the" of a matrix. This value
    #       is computed as the square root of the maximum eigenvalue
    #       of A^H A where A^H is the conjugate transpose.
    #   For p = Inf:
    #       The maximum absolute row sum norm |A|_inf is defined
    #       as the maximum of the sum of the absolute valued elements
    #       of rows of the matrix.

    # FUNCTION:

    # Compute Norm:
    ans <- NA
    if (p == 1) {
        x <- abs(x)
        ans <- max(apply(x, 2, sum))
    }
    if (p == 2) {
        ans <- sqrt(max(eigen(t(x) %*% x)$values))
    }
    if (p == Inf) {
        x <- abs(x)
        ans <- max(apply(x, 1, sum))
    }
    if (is.na(ans)) stop("Invalid value for p")

    # Return value:
    ans
}


################################################################################

