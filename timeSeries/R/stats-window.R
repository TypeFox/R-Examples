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
#  window,timeSeries         Extracts a piece from a 'timeSeries' object
# DEPRECATED:               DESCRIPTION:
#  cut,timeSeries            Extracsts a piece from a 'timeSeries' object
################################################################################


.window.timeSeries <- 
    function(x, start, end, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Windows a piece from a 'timeSeries' object.

    # Arguments:
    #   x - a 'timeSeries' object
    #   from, to - two 'timeDate' position vectors which size the
    #       blocks

    # Details:
    #   from and to, are both included in the window.

    # Value:
    #   Returns a S4 object of class 'timeSeries'.

    # FUNCTION:
      
    # Check Arguments:
    stopifnot(is.timeSeries(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
      
    # Check for Signal Series
    if (x@format == "counts")
        stop(as.character(match.call())[1],
             " is for time series and not for signal series.")

    # check if all argument names are used
    if (length(dot <- list(...))) {
        if (any(names(dot) %in% c("from", "to"))) {
            if (!is.null(from <- dot$from)) start <- from
            if (!is.null(to <- dot$to)) end <- to
            warning("Arguments 'from/to' are deprecated.\nUse instead 'start/end'.", call. = FALSE)
        }
    }

    start <- timeDate(start)
    end <- timeDate(end)
    Positions <- time(x)
    test <- (Positions >= start & Positions <= end)
    ans <- x[test,]
      
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation 
      
    # Return Value:
    ans
}


setMethod("window", "timeSeries",
          function(x, start, end, ...) .window.timeSeries(x, start, end, ...))

          
# until UseMethod dispatches S4 methods in 'base' functions
window.timeSeries <- function(x, ...) .window.timeSeries(x, ...)


###############################################################################


.cut.timeSeries <- 
    function (x, from, to, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Cuts out a piece from a 'timeSeries' object.

    # Arguments:
    #   x - a 'timeSeries' object
    #   from, to - two 'timeDate' position vectors which size the
    #       blocks

    # Value:
    #   Returns a S4 object of class 'timeSeries'.

    # FUNCTION:

    # .Deprecated("window", "timeSeries")

    stopifnot(is.timeSeries(x))
    if (x@format == "counts")
        stop(as.character(match.call())[1],
             " is for time series and not for signal series.")

    from = timeDate(from)
    to = timeDate(to)
    Positions = time(x)

    test = (Positions >= from & Positions <= to)
    ans <- x[test,]

    # Return value:
    ans
}


setMethod("cut", "timeSeries",
    function (x, from, to, ...) .cut.timeSeries(x, from, to, ...))

          
# until UseMethod dispatches S4 methods in 'base' functions
cut.timeSeries <- function(x, ...) .cut.timeSeries(x, ...)


################################################################################

