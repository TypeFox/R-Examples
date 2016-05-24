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
#  dummySeries               Creates a dummy monthly 'timeSeries' object
#  dummyDailySeries          Creates a dummy daily 'timeSeries' object
################################################################################


# DW:
# A more natural name for the function dummySeries() would be 
# dummyMonthlySeries() to have the same naming conventions like in the
# case of the dummy daily series.


dummySeries <- 
  function(...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a monthly dummy 'time Series' object
    
    # Arguments:
    #   ... - optional arguments passed to the function timeSeries().
    
    # FUnction:
    
    # Return Value:
    timeSeries(matrix(runif(24), ncol = 2), as.character(timeCalendar()), ...)
}

    
# ------------------------------------------------------------------------------


dummyDailySeries <-
  function(x = rnorm(365), units = NULL, zone = "", FinCenter = "")
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a dummy daily time Series

    # Arguments:
    #   x - a numeric vector
    #   origin - the first date in the series

    # FUNCTION:
    if (zone == "")
        zone <- getRmetricsOptions("myFinCenter")
    if (FinCenter == "")
        FinCenter <- getRmetricsOptions("myFinCenter")

    # Check:
    stopifnot(is.numeric(x))
    if (is.null(units)) units <- paste("X", 1:NCOL(x), sep = "")
    stopifnot(length(units) == NCOL(x))

    # Time Series:
    if (is.vector(x)) data = matrix(x, ncol = 1)
    if (is.matrix(x)) data = x
    positions <- timeSequence(from = "1970-01-01", length.out =
        NROW(data), zone = zone, FinCenter = FinCenter)
    ans <- timeSeries(data = data, charvec = positions, units = units,
        zone = zone, FinCenter = FinCenter)

    # Return Value:
    ans
}


################################################################################

