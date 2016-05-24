
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
#  julian,timeDate           Returns Julian day counts since 1970-01-01
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("julian", "timeDate",
    function(x, origin = timeDate("1970-01-01"),
    units = c("auto", "secs", "mins", "hours", "days", "weeks"),
    zone = NULL, FinCenter = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Extracts Julian time in days since 1970-01-01

    # Arguments:
    #   x - a 'timeDate' object
    #   units - a character string, one of the units listed,
    #       by default "secs".

    # Value:
    #   Returns the number of days (possibly fractional) since
    #   the origin.

    # Details:
    #   The origin is "1970-01-01 00:00:00 GMT"

    # Set Timezone to GMT:

    # Check Class Type:
    stopifnot(is(x, "timeDate"))
    units = match.arg(units)

    # POSIX:
    if (is.null(zone)) zone = x@FinCenter
    if (is.null(FinCenter)) FinCenter = x@FinCenter
    ct = timeDate(x, zone = zone, FinCenter = FinCenter)

    # Difftime:
    if (is.null(origin))
        origin = timeDate("1970-01-01", zone = "GMT", FinCenter = "GMT")
    res = difftimeDate(ct, origin, units = units)

    # Return Value:
    structure(res, origin = origin)
})


################################################################################

