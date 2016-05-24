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


###############################################################################
# FUNCTION:                 DESCRIPTION:
#  cumulated                 Computes cumulated series from financial returns
#  cumulated.default         Computes cumulated series, default method
###############################################################################


cumulated <-
  function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes cumulated series from financial returns

    # Return Value:
    UseMethod("cumulated")
}


# ------------------------------------------------------------------------------


cumulated.default <-
  function(x, method = c("continuous", "discrete", "compound", "simple"),
    percentage = FALSE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes cumulated series from financial returns
    #   supports 'matrix' and 'timeSeries'.

    # Arguments:
    #   x - data object containing ordered price observations
    #   method - "continuous == "compound" and "discrete" == "simple"

    # Example:
    #   X = as.timeSeries(data(msft.dat))[1:10, "Close"]; X = X/series(X)[1, 1]
    #   x = returns(X, "continuous"); x;  X; cumulated(x, "continuous")
    #   x = returns(X, "discrete"); x;  X; cumulated(x, "discrete")

    # Note:
    #   To make it conform with PortfolioAnalytics:
    #   "compound" == "continuous", and "simple" == "discrete"

    # FUNCTION:

    # Check Arguments:
    stopifnot(is.timeSeries(x))
    
    # Settings:
    method <- match.arg(method)
    
    # Extract Title and Documentation:
    Title <- x@title
    Documentation <- x@documentation
    
    # Handle Missing Values:
    # if (na.rm) x = na.omit(x, ...)

    # Transform data:
    if (percentage) x <- x/100
    positions <- time(x)

    # Calculate Cumulates:
    # ... colCumsums and colCumprods are generic functions with
    #     methods for 'matrix' and 'timeSeries'.
    if(method == "geometric") {
        ans <- colCumsums(x)
    }
    if(method == "compound" || method == "continuous") {
        ans <- exp(colCumsums(x))
    }
    if(method == "simple" || method == "discrete") {
        ans <- colCumprods(1+x)
    }
    
    # Preserve Title and Documentation:
    ans@title <- Title
    ans@documentation <- Documentation

    # Return Value:
    ans
}


################################################################################

