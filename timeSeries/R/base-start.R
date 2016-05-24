#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# FUNCTION:                 DESCRIPTION:
#  start,timeSeries          Extracts start date of a 'timeSeries' object
#  end,timeSeries            Extracts end date of a 'timeSeries' object
################################################################################


.start.timeSeries <- function(x, ...)
{
    # Description:
    #   Extracts start date of a 'timeSeries' object
    
    # FUNCTION:
    
    # Extract Date:
    if (length(x@positions)>0)
        timeDate(min(x@positions), zone = "GMT", FinCenter = x@FinCenter)
    else
        NULL
}

setMethod("start" , "timeSeries", function(x, ...) .start.timeSeries(x, ...))

# until UseMethod dispatches S4 methods in 'base' functions
start.timeSeries <- function(x, ...) .start.timeSeries(x, ...)


# ------------------------------------------------------------------------------


.end.timeSeries <- function(x, ...)
{
    # Description:
    #   Extracts start date of a 'timeSeries' object
    
    # FUNCTION:
    
    # Extract Date:
    if (length(x@positions)>0)
        timeDate(max(x@positions), zone = "GMT", FinCenter = x@FinCenter)
    else
        NULL
}

setMethod("end", "timeSeries", function(x, ...) .end.timeSeries(x, ...))

# until UseMethod dispatches S4 methods in 'base' functions
end.timeSeries <- function(x, ...) .end.timeSeries(x, ...)


################################################################################

