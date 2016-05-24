
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
# MEHODS:                   DESCRIPTION:
#  diff.timeDate             Returns suitably lagged and iterated differences
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
diff.timeDate <-
    function (x, lag = 1, differences = 1, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns suitably lagged and iterated differences

    # Arguments:
    #   x - a 'timeDate' object.
    #   lag - an integer indicating which lag to use, by
    #       default 1.
    #   differences - an integer indicating the order of the
    #       difference, by default 1.
    #   ... - further arguments to be passed to or from methods.

    # Value:
    #   If 'x' is a vector of length 'n' and 'differences=1', then
    #   the computed result is equal to the successive differences
    #   'x[(1+lag):n] - x[1:(n-lag)]'.
    #   If 'difference' is larger than one this algorithm is applied
    #   recursively to 'x'. Note that the returned value is a vector
    #   which is shorter than 'x'.

    # FUNCTION:

    # Convert to GMT:
    GMT <- timeDate(x, zone = x@FinCenter, FinCenter = "GMT")

    diff(as.POSIXct(GMT@Data),
         lag = lag, differences = differences, ...)
}


################################################################################

