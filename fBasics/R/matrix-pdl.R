
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
#  pdl                       Regressor matrix for polynomial distributed lags
################################################################################


pdl <-
function(x, d = 2, q = 3, trim = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Regressor matrix for polynomial distributed lags

    # Aruments:
    #   x - a numeric vector.
    #   d - an integer specifying the order of the polynomial.
    #   q - an integer specifying the number of lags to use in
    #       creating polynomial distributed lags. This must be
    #       greater than d.
    #   trim - a logical flag; if TRUE, the missing values at
    #       the beginning of the returned matrix will be trimmed.

    # Value:
    #   Returns a matrix representing the regressor matrix.

    # Example:
    #   stack.loss = c(
    #       42, 37, 37, 28, 18, 18, 19, 20, 15, 14, 14,
    #       13, 11, 12,  8,  7,  8,  8,  9, 15, 15)
    #   pdl(stack.loss)

    # FUNCTION:

    # Check:
    stopifnot(q > d)

    # Polynomial distributed lags:
    M = tslag(x, 1:q, FALSE)
    C = NULL
    for (i in 0:d) { C = rbind(C, (1:q)^i) }
    Z = NULL
    for (i in 1:(d+1)) { Z = cbind(Z, apply(t(C[i,]*t(M)), 1, sum)) }
    Z[, 1] = Z[, 1] + x

    # Trim:
    if (trim) {
        indexes = (1:length(Z[,1]))[!is.na(apply(Z, 1, sum))]
        Z = Z[indexes, ] }

    # Return Value:
    Z
}


################################################################################

