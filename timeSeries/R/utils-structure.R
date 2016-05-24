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
# FUNCTION:               DESCRIPTION:
#  str,timeSeries          Displays the structure of a 'timeSeries' object
################################################################################


.str.timeSeries <- 
    function(object, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Displays the structure of a 'timeSeries' object.

    # Arguments:
    #   object - an object of class 'timeSeries'.
    #   ... - 

    # FUNCTION:

    # Series Name:
    cat("Time Series:          ")
    cat("\n Name:              ", as.character(c(substitute(object))))
    
    # YC : as.character(c( important to handle str(timeSeries())
    
    # Data Matrix:
    Dim = dim(object)
    cat("\nData Matrix:        ")
    cat("\n Dimension:         ", Dim)
    cat("\n Column Names:      ", colnames(object) )
    firstName = rownames(object)[1]
    lastName = rownames(object)[Dim[1]]
    cat("\n Row Names:         ", firstName, " ... ", lastName)
    
    # Date/Time Positions:
    cat("\nPositions:          ")
    cat("\n Start:             ", as.character(start(object)))
    cat("\n End:               ", as.character(end(object)))
    
    # Other Attributes:
    cat("\nWith:               ")
    cat("\n Format:            ", object@format)
    cat("\n FinCenter:         ", object@FinCenter)
    cat("\n Units:             ", object@units)
    cat("\n Title:             ", object@title)
    cat("\n Documentation:     ", object@documentation)
    cat("\n")

    # Return Value:
    invisible()
}


setMethod("str", "timeSeries",
    function(object, ...) .str.timeSeries(object, ...))

    
# until UseMethod dispatches S4 methods in 'base' functions
str.timeSeries <- function (object, ...) .str.timeSeries(object, ...)


################################################################################

