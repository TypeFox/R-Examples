
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
#  window.timeDate           Extracts a piece from a 'timeDate' object
#  cut.timeDate              Extracts a piece from a 'timeDate' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
window.timeDate <-
    function(x, start, end, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Descriptions:
    #   Subsets a timeDate object between from and to dates.

    # Arguments
    #   x - an object of class "timeDate"
    #   start, end - start and end dates as "timeDate" objects.

    # Note:
    #   This is synonome for the function window() which should
    #   be preferred.

    # FUNCTION:

    # Extract Subset:
    ct <- as.POSIXct(x)
    FROM = as.POSIXct(timeDate(start, zone = x@FinCenter, FinCenter = "GMT"))
    TO = as.POSIXct(timeDate(end, zone = x@FinCenter, FinCenter = "GMT"))
    test = (ct >= FROM & ct <= TO)

    if (all(!test)) {
        ans <- x
        ans@Data <- as.POSIXct(NA)
    } else {
        ans <- x[test]
    }

    # Return Value:
    ans
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
cut.timeDate <-
    function(x, from , to, ...)
{
    # A function implemented by Diethelm Wuertz

    # Descriptions:
    #   Subsets a timeDate object between from and to dates.

    # Arguments
    #   x - an object of class "timeDate"
    #   from, to - start and end dates as "timeDate" objects.

    # Note:
    #   This is synonome for the function window() which should
    #   be preferred.

    # FUNCTION:

    # Extract Subset:
    X = timeDate(x, zone = x@FinCenter, FinCenter = "GMT")
    FROM = timeDate(from, zone = x@FinCenter, FinCenter = "GMT")
    TO = timeDate(to, zone = x@FinCenter, FinCenter = "GMT")
    test = (X >= FROM & X <= TO)
    ans = timeDate(X[test], zone = "GMT", FinCenter = x@FinCenter)

    # Return Value:
    ans
}


################################################################################

