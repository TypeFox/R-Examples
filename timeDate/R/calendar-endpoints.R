
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
# FUNCTION:              DESCRIPTION:
#  endpoints              Returns endpoint indexes from a timeDate object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
## YC: do not make this function visible unless one rename it to avoid
## conflicts with xts endpoints() function.
.endpoints <-
function (x, on = c("months", "years", "quarters", "weeks", "days",
    "hours", "minutes", "seconds"), k = 1)
{
    # Description:
    #   Returns endpoint indexes from a 'timeDate' object.

    # Arguments:
    #   x - a 'timeDate' object
    #   on - the periods endpoints to find as a character string
    #   k - along every k-th element

    # Note:
    #   Behaves like function entpoints() from R package xts, which
    #     extracts index values of a given "xts" object corresponding
    #     to the last observations given a period specified by "on",
    #   with a zero added to the beginning of the vector, and the
    #     index of the last observation in x at the end.
    #   The Index rules are borrowed from Jeff Ryans endpoints()
    #     function.

    # FUNCTION:

    # Convert to POSIX:
    on <- match.arg(on)
    posix <- as.POSIXct(x)
    .posix <- unclass(posix)

    # Apply index extraction rules:
    if (on == "years") {
        ans <- as.integer(which(diff(as.POSIXlt(posix)$year %/% k + 1) != 0))
    } else if (on == "quarters") {
        ans <- as.integer(which(diff((as.POSIXlt(posix)$mon %/% 3) + 1) != 0))
    } else if (on == "months") {
        ans <- as.integer(which(diff(as.POSIXlt(posix)$mon %/% k + 1) != 0))
    } else if (on == "weeks") {
        ans <- as.integer(which(diff((.posix + (3L * 86400L))%/%604800L %/% k + 1) != 0))
    } else if (on == "days") {
        ans <- as.integer(which(diff(.posix %/% 86400L %/% k + 1) != 0))
    } else if (on == "hours") {
        ans <- as.integer(which(diff(.posix %/% 3600L %/% k + 1) != 0))
    } else if (on == "minutes" || on == "mins") {
        ans <- as.integer(which(diff(.posix %/% 60L %/% k + 1) != 0))
    } else if (on == "seconds" || on == "secs") {
        ans <- as.integer(which(diff(.posix %/% k + 1) != 0))
    }
    ans <- c(0, ans, NROW(x))

    # Return Value:
    ans
}


################################################################################
