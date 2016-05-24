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
# METHOD:                   SUBSETTING METHODS ON DATA:
#  head,timeSeries           Returns the head of a 'timeSeries' object
#  tail,timeSeries           Returns the tail of a 'timeSeries' object
################################################################################


.head.timeSeries <- 
    function(x, n = 6, recordIDs = FALSE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the head of a 'timeSeries' object

    # Arguments:
    #   x - a 'timeSeries' object.
    #   n - a single integer. If positive, number of the first n records (rows)
    #        to  be returned. If negative, all but the n first number of 
    #        elements of x are returned.
    #   recordIDs - a logical flag, should the record identification
    #       be shown? By default FALSE.
    #   ... - 

    # Value:
    #   Returns the tail of an object of class 'timeSeries'.
    
    # FUNCTION:

    # Head:
    if (recordIDs & dim(x)[1] == dim(x@recordIDs)[1])
        cbind(head.matrix(x, n = n, ...), head(x@recordIDs, n = n, ...))
    else
        head.matrix(x, n = n, ...)
}


setMethod("head", "timeSeries", function(x, n = 6, recordIDs = FALSE, ...)
          .head.timeSeries(x, n, recordIDs, ...))

          
# until UseMethod dispatches S4 methods in 'base' functions
head.timeSeries <- function(x, ...) .head.timeSeries(x, ...)


# ------------------------------------------------------------------------------


.tail.timeSeries <- 
    function(x, n = 6, recordIDs = FALSE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the tail of a 'timeSeries' object

    # Arguments:
    #   x - a 'timeSeries' object.
    #   n - a single integer. If positive, number of the last n records (rows)
    #       to be returned. If negative, all but the n last number of 
    #       elements of x are returned.    
    #   recordIDs - a logical flag, should the record identification
    #       be shown? By default FALSE.
    #   ... - 

    # Value:
    #   Returns the tail of an object of class 'timeSeries'.

    # FUNCTION:

    # Tail:
    if (recordIDs & dim(x)[1] == dim(x@recordIDs)[1])
        cbind(tail.matrix(x, n = n, addrownums = FALSE, ...),
              tail(x@recordIDs, n = n, addrownums = FALSE, ...))
    else
        tail.matrix(x, n = n, addrownums = FALSE, ...)
}


setMethod("tail", "timeSeries", function(x, n = 6, recordIDs = FALSE, ...)
          .tail.timeSeries(x, n, recordIDs, ...))

          
# until UseMethod dispatches S4 methods in 'base' functions
tail.timeSeries <- function(x, ...) .tail.timeSeries(x, ...)


################################################################################

