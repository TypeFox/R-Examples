
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
# Function:               DESCRIPTION:
#  vec                     Stacks a matrix as column vector
#  vech                    Stacks a lower triangle matrix
################################################################################


vec <-
function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Stacks a matrix as column vector

    # Details:
    #   vec(X) = (X11, X21, ..., XN1, X12, X22, ..., XNN)'

    # Note:
    #   Example for a 3x3 Matrix:
    #   X11, X21, X22, X31, X32, X33

    # FUNCTION:

    # Return Value:
    t(t(as.vector(x)))
}


# ------------------------------------------------------------------------------


vech <-
function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Stacks a lower triangle matrix

    # Details:
    #   vech is the operator that stacks the lower triangle
    #   of a NxN matrix as an N(N+1)/2x1 vector:
    #   vech(X) =(X11, X21, X22, X31, ..., XNN)'

    # Note:
    #   Example for a 3x3 Matrix:
    #   X11, X21, X22, X31, X32, X33

    # FUNCTION:

    # Return Value:
    t(x[!upper.tri(x)])
}


################################################################################

