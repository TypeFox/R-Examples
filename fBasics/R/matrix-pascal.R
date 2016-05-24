
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
#  pascal                    Creates a Pascal matrix
################################################################################


pascal <-
function(n)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a Pascal matrix

    # Arguments:
    #   n - the dimension of the square matrix

    # Details:
    #   http://mathworld.wolfram.com/PascalMatrix.html
    #   Pascal matrices are symmetric and positive definite.
    #   The determinant of a Pascal matrix is 1.
    #   The inverse of a Pascal matrix has integer entries.
    #   If lambda is an eigenvalue of a Pascal matrix,
    #       then 1/lambda is also an eigenvalue of the matrix.
    #   The Cholesky factor of a Pascal matrix consists of
    #       the elements of Pascal's triangle

    # FUNCTION:

    # Pascal:
    N = n-1
    n.over.r = function(n, r) {
        prod(1:n) / (prod(1:(n-r)) * prod(1:r) ) }
    X = rep(1, N)
    for ( i in 1:N )
        for ( j in 1:N )
        X = c(X, n.over.r(i+j, j))
        X = cbind(rep(1, N+1), matrix(X, byrow = TRUE, ncol = N))

    # Return Value:
    X
}


################################################################################

