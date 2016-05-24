
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

# Copyrights (C)
# for this R-port:
#   1999 - Diethelm Wuertz, GPL
#   2007 - Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@phys.ethz.ch>
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# MEHODS:                   DESCRIPTION:
#  rep.timeDate              Replicates a 'timeDate' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
rep.timeDate <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Replicates objects of class 'timeDate'

    # Arguments:
    #   x - a 'timeDate' object
    #   times - a non-negative integer.  A vector giving the number
    #       of times to repeat each element if of length 'length(x)',
    #       or to repeat the whole vector if of length 1.

    # Value:
    #   Returns a vector of repeated elements belonging to the same
    #   class as 'x'.

    # FUNCTION:

    # Replicate:
    ct <- rep(as.POSIXct(x), ...)
    ans <- timeDate(ct, zone = "GMT", FinCenter = x@FinCenter)

    # Return Value:
    ans
}


################################################################################

