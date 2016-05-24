
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
#  tr                        Returns the trace of a matrix
################################################################################


tr <-
function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the trace of a matrix

    # FUNCTION:

    # Trace:
    if (dim(x)[1] != dim(x)[2] ) {
        return(NA)
    } else {
        return(sum(diag(x)))
    }

    # Return Value:
    invisible()
}


################################################################################

