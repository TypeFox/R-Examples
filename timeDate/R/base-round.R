
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
#  round.timeDate            Rounds objects of class 'timeDate'
#  trunc.timeDate            Truncates objects of class 'timeDate'
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
round.timeDate <- function(x, digits = c("days", "hours", "mins"))
{
    # A function implemented by Diethelm Wuertz
    # and modified by Yohan Chalabi

    # Note:
    #   round.timeDate(x, units = c("days", "hours", "mins"))    # FAILS !!!

    # FUNCTION:

    # Get Units:
    units <- match.arg(digits)
    FinCenter <- finCenter(x)

    # Use:
    lt <- round.POSIXt(as.POSIXlt(x, tz = "GMT"), units = units)
    ans <- timeDate(lt, zone = FinCenter, FinCenter = FinCenter)

    # Return Value:
    ans
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
trunc.timeDate <- function(x, units = c("days", "hours", "mins"), ...)
{
    # A function implemented by Diethelm Wuertz
    # and modified by Yohan Chalabi

    # FUNCTION:

    # Get Units:
    units = match.arg(units)
    FinCenter <- finCenter(x)

    # Sorting under GMT is not what we want!
    # GMT = timeDate(x, zone = x@FinCenter, FinCenter = "GMT")
    # lt = trunc.POSIXt(GMT@Data, units = units)
    # ans = timeDate(lt, zone = "GMT", FinCenter = x@FinCenter)

    # Use:
    lt <- trunc.POSIXt(as.POSIXlt(x, tz = "GMT"), units = units)
    ans <- timeDate(lt, zone = FinCenter, FinCenter = FinCenter)

    # Return Value:
    ans
}


################################################################################

