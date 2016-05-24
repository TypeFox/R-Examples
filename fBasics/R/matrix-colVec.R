
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
#  colVec                    Creates a column vector from a data vector
#  rowVec                    Creates a row vector from a data vector
################################################################################


colVec <-
function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Converts a vector to a column vector

    # Details:
    #   A column vector is a matrix with one column.

    # Return Value:

    # FUNCTION:

    # Double Transpose:
    ans = t(t(x))

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


rowVec <-
function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Converts a vector to a row vector

    # Details:
    #   A row vector is a matrix with one row.

    # FUNCTION:

    # Transpose:
    ans = t(x)

    # Return Value:
    ans
}


################################################################################

