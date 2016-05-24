
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
#  runlengths                Returns 'timeSeries' object of runlengths
################################################################################


runlengths <-
    function(x, ...)
{
    # A function implemetned by Diethelm Wuertz

    # Description:
    #   Returns 'timeSeries' object of runlengths

    # Arguments:
    #   x - an univariate 'timeSeries' object of financial returns
    #   ... - arguments passed to the function na.omit()

    # Value:
    #   runlengths an object of class 'timeSeries'.

    # Note:
    #   Zeroes are handled as NA.

    # Example:
    #   set.seed(4711)
    #   x.tS = timeSeries(data=rnorm(12), charvec=timeCalendar(), units="x")
    #   runlengths(x.tS)

    # FUNCTION:

    # Check arguments:
    stopifnot(is.timeSeries(x))
    stopifnot(isUnivariate(x))
      
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
      
    # Handle Missing Values:
    x[x == 0] <- NA
    x.vec = sign(as.vector(na.omit(x, ...)))

    # Compute Run Lengths:
    n <- length(x.vec)
    y <- x.vec[-1L] != x.vec[-n]
    Index <- c(which(y | is.na(y)), n)
    X = x[Index, ]
    series(X) <- matrix(diff(c(0L, Index)), ncol = 1)

    # Reset recordIDs:
    X@recordIDs <- data.frame()
    
    # Preserve Title and Documentation:
    X@title <- Title
    X@documentation <- Documentation
      
    # Return Value:
    X
}


################################################################################

