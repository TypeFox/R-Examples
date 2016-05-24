
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
# FUNCTION:                      DESCRIPTION:
#  isDaily,timeDate-method        Tests 'timeDate' has daily time stamps
#  isMonthly,timeDate-method      Tests 'timeDate' has monthly time stamps
#  isQuarterly,timeDate-method    Tests 'timeDate' has quarterly time stamps
#  isRegular,timeDate-method      Tests 'timeDate' has regular time stamps
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("isDaily", "timeDate", function(x)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Test if a timeDate object has daily time stamps
    
    # Example:
    #   isDaily(timeSequence(by = "day", length.out = 20))
    #   isDaily(timeCalendar())
    #   isDaily(timeSequence(by = "hour", length.out = 100))
    
    # Details:
    #   Definition: A timeDate Object is a Daily timeDate object
    #   if we have not more than one date/time stamp per day.
    
    # Arguments:
    #   x - an object of class timeDate
    
    # FUNCTION:
    
    # Daily ?
    num <- as.numeric(as.POSIXct(x))
    daily <- seq(from = num[1], by = 60*60*24, length.out=length(num))
    
    # Return
    (identical(daily, num))
})


################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("isMonthly", "timeDate", function(x)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Tests if a timeDate object has monthly time stamps
    
    # Arguments:
    #   x - an object of class timeDate
    
    # Details:
    #   Definition: A timeDate Object is a Monthly timeDate object
    #   if we have not more than one date/time stamp per month.
    #   Note a monthly series is also a daily series.
    
    # Example:
    #   isMonthly(timeSequence(by = "day", length.out = 20))
    #   isMonthly(timeCalendar())
    #   isDaily(timeCalendar())
    #   isMonthly(timeSequence(by = "hour", length.out = 100))
    
    # FUNCTION:
    
    # Monthly ?
    m <- c(timeDate::months(x)) #-> c() to remove attributes
    # (m[1] -1) -> shift vector to match first entry in m
    monthly <- seq(from = m[1]-1, length.out=length(m)) %% 12 + 1
    
    # Return
    (identical(monthly, m))
})


################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("isQuarterly", "timeDate", function(x)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Tests if a timeDate object has quarterly time stamps
    
    # Arguments:
    #   x - an object of class timeDate
    
    # Details:
    #   Definition: A timeDate Object is a Quarterly timeDate object
    #   if we have not more than one date/time stamp per quarter
    #   Note a quarterly series is also a daily and a monthly series.
    
    # Example:
    #   isQuarterly(timeSequence(by = "day", length.out = 20))
    #   isQuarterly(timeCalendar())
    #   isQuarterly(timeSequence(by = "hour", length.out = 100))
    #   isQuarterly(timeCalendar()[(1:4)*3])
    #   isMonthly(timeCalendar()[(1:4)*3])
    #   isDaily(timeCalendar()[(1:4)*3])
    
    # FUNCTION:
    
    # Quartertly ?
    m <- c(timeDate::months(x)) #-> c() to remove attributes
    # (m[1] -1) -> shift vector to match first entry in m
    quarterly <- seq(from = m[1]-1, by = 3, length=length(m)) %% 12 + 1
    
    # Return
    (identical(quarterly, m))
})


################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("isRegular", "timeDate",  function(x)
{
    # A function implemented by Diethelm Wuertz
    
    # Descriptions:
    #   Tests if a timeDate object has regular time stamps
    
    # Example:
    #   isRegular(timeSequence(by = "day", length.out = 20))
    #   isRegular(timeCalendar())
    #   isRegular(timeSequence(by = "hour", length.out = 100))
    
    # Details:
    #   Definition: A timeDate Object is a Regular timeDate object
    #   if the timeDate object is either monthly or quarterly,
    #   otherwise not.
    
    # Arguments:
    #   x - an object of class timeDate
    
    # FUNCTION:
    
    # Regular ?
    (isMonthly(x) | isQuarterly(x))
    })


################################################################################