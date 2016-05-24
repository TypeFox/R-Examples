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
#  sort,timeSeries           Sorts a 'timeSeries' object in time
################################################################################


.sort.timeSeries <- function (x, decreasing = FALSE, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Time sorts a 'timeSeries' object

    # Arguments:
    #   x - a 'timeSeries' object.

    # Value:
    #   Returns a time sorted object of class 'timeSeries'.

    # FUNCTION:

    # check if really necessary to sort x
    # important in order to improve efficiency
    ## NB: is.unsorted can return NA
    if (!decreasing && !isTRUE(is.unsorted(x@positions))) return(x)

    if (length(x@positions)>0)
        x[order(x@positions, decreasing = decreasing), ]
    else
        x
}

setMethod("sort", "timeSeries", function (x, decreasing = FALSE, ...)
          .sort.timeSeries(x, decreasing = decreasing, ...))

# until UseMethod dispatches S4 methods in 'base' functions
sort.timeSeries <- function(x, decreasing = FALSE, ...)
    .sort.timeSeries(x, decreasing = decreasing, ...)

################################################################################

