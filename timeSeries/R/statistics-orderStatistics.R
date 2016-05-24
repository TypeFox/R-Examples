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
#  orderStatistics           Compute order statistic of a 'timeSeries' object
################################################################################


orderStatistics <-
function(x)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Compute the order statistics for a 'timeSeries object

    # Value:
    #   A named list with the order statistics for each column of
    #   the inputted series.

    # FUNCTION:

    # Order Statistics:
    td <- time(x)
    
    # Return Value:
    mapply(
        function(cl, nm) {
           S <- sort(cl, index.return = TRUE)
           timeSeries(data = S$x, charvec = td[S$ix], units = nm)}, 
        as.list(x), 
        colnames(x), 
        SIMPLIFY = FALSE)
}


################################################################################

