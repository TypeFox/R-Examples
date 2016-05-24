
# This R package is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This R package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this R package; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# METHOD:                   DESCRIPTION:
#  atoms,timeDate            Returns date/time atoms from a 'timeDate' object
#  months,timeDate           Extracts months atom from a 'timeDate' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("atoms", "timeDate", function(x, ...)
{
    # A function implemented by Diethelm Wuertz
    # and improved by Yohan Chalabi

    # Description:
    #   Extracts atoms from a 'timeDate' object.

    # Arguments:
    #   x - a 'timeDate' object from which to extract the
    #       calendar "atoms".

    # Value:
    #   Returns a data.frame with the following calendar atoms:
    #   Y(ear), m(onth), d(ay), H(our), M(inutes), S(econds).

    # FUNCTION:

    # Check Class Type:
    if (!inherits(x, "timeDate")) stop("Wrong class type")

    # mdy:
    X <- as.POSIXlt(x, tz = "GMT")
    Y <- X$year + 1900
    m <- X$mon + 1
    d <- X$mday
    H <- X$hour
    M <- X$min
    S <- X$sec

    # Data Frame:
    ans <- data.frame(Y = Y, m = m, d = d, H = H, M = M, S = S)
    attr(ans, "control") <- c(FinCenter = finCenter(x))

    # Return Value:
    ans
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("atoms", "ANY", function(x, ...)
{
    # A function implemented by Diethelm WUertz

    # FUNCTION:

    # Return Value:
    invisible(x)
})


################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("months", "timeDate",
    function(x, abbreviate=NULL)
{
    # A function implemented by Diethelm Wuertz
    # and improved by Yohan Chalabi

    # Description:
    #   Extracts months atom from a timeDate object

    # Arguments:
    #   x - a 'timeDate' object from which to extract the
    #       month "atom".

    # Value:
    #   Returns the month from a 'timeDate' object as an integer
    #   value or vector with elements ranging between 1 and 12,
    #   numbering the months from January to December.

    # FUNCTION:

    # Check Class Type:
    if (!inherits(x, "timeDate")) stop("Wrong class type")

    # Month:
    ans <- as.POSIXlt(x, tz = "GMT")$mon+1
    attr(ans, "control") <- c(FinCenter = finCenter(x))

    # Return Value:
    ans
})


################################################################################

