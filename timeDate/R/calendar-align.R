
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
# FUNCTION:               DESCRIPTION:
#  align,timeDate          Aligns a 'timeDate' object to regular time stamps
#  align,ANY
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("align", "timeDate",
    function(x, by = "1d", offset = "0s")
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi
    
    # Description:
    #   Aligns a 'timeDate' object to regular time stamps

    # Example:
    #   align(timeCalendar(), "1w")        # Weekly
    #   align(timeCalendar(), "2w", "3d")  # Bi-Weekly with offset

    # FUNCTION:

    # Settings:
    periods <- c(7*24*3600, 24*3600, 3600, 60, 1)
    names(periods) <- c("w", "d", "h", "m", "s")
    offset <- as.integer(gsub("[a-z]", "", offset, perl = TRUE)) *
        periods[gsub("[ 0-9]", "", offset, perl = TRUE)]
    offset <- as.vector(offset)

    # Return Value:
    seq(from = x[1] + offset, to = x[length(x)], by = by)

})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("align", "ANY",
    function(x, y, xout, method = "linear", n = 50, rule = 1, f = 0,
    ties = mean, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # FUNCTION:

    # Align by Approximation:
    ans = approx(x = x, y = y, xout = xout, method = method, n = n,
        rule = rule, f = f, ties = ties, ...)

    # Return Value:
    ans
})


################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
alignDaily <-
function(x, include.weekends=FALSE)
{
     # A function implemented by Diethelm Wuertz
     
     # Description:
     #    Aligns a 'timeDate' object to end-of-day dates
     
     # Arguments:
     #    x - a 'timeDate' object
     #    include.weekends - a logical, should weekends be included?
     
     # FUNCTION:
     
     # Align:
     if (include.weekends) {
         tD <- align(x)
     } else {
         tD <- align(x)
         tD <- tD[isWeekday(tD)]
     }
     
     # Return Value:
     tD
}


# ----------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
alignMonthly <- 
function(x, include.weekends=FALSE)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Aligns a 'timeDate' object to end of month dates
     
    # Arguments:
    #    x - a 'timeDate' object
    #    include.weekends - a logical, should weekends be included?
     
    # FUNCTION:
    
    # Align:
    if (include.weekends) {
        tD <- timeLastDayInMonth(x)
    } else {
        tD <- timeLastDayInMonth(x)
        tD[isWeekend(tD)] <- tD[isWeekend(tD)] - 24*3600
    }
    
    # Return Value:
    tD
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
alignQuarterly <- 
function(x, include.weekends=FALSE)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #    Aligns a 'timeDate' object to end-of-quarter dates
    
    # Arguments:
    #    x - a 'timeDate' object
    #    include.weekends - a logical, should weekends be included?
     
    # FUNCTION:
    
    # Align:
    if (include.weekends) {
        tD <- timeLastDayInQuarter(x)
    } else {
        tD <- timeLastDayInQuarter(x)
        tD[isWeekend(tD)] <- tD[isWeekend(tD)] - 24*3600
    }
    
    # Return Value:
    tD
}


###############################################################################

