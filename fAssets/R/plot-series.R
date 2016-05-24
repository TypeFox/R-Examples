
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                   DESCRIPTION:
#  assetsReturnPlot            Displays time series of individual assets
#  assetsCumulatedPlot         Displays time series of individual assets
#  assetsSeriesPlot            Displays time series of individual assets
################################################################################


assetsReturnPlot =
    function(x, col = "steelblue", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays return series of individual assets

    # Arguments:
    #   x - a timeSeries object of financial returns or any other 
    #       rectangular object which can be transformed by the 
    #       function as.matrix into a numeric matrix.

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))
    #   par(mfrow = c(3,3)); assetsReturnPlot(x); par(mfrow = c(1,1))
    
    # FUNCTION:

    # Settings:
    n = ncol(x)
    if (length(col) == 1) col = rep(col, times = n)

    # Plot:
    seriesPlot(x, ylab = "Returns", col = col, ...)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


assetsCumulatedPlot =
    function(x, col = "steelblue", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays cumulated return series of individual assets

    # Arguments:
    #   x - a timeSeries object of financial returns or any other 
    #       rectangular object which can be transformed by the 
    #       function as.matrix into a numeric matrix.

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))
    #   par(mfrow = c(3,3)); assetsCumulatedPlot(x); par(mfrow = c(1,1))
    
    # FUNCTION:

    # Settings:
    n = ncol(x)
    if (length(col) == 1) col = rep(col, times = n)

    # Plot:
    x = exp(colCumsums(x))
    seriesPlot(x, ylab = "Cumulated Returns", col = col, ...)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


assetsSeriesPlot =
    function(x, col = "steelblue", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays a derived series of individual assets

    # Arguments:
    #   x - a timeSeries object or any other rectangular object
    #       which can be transformed by the function as. matrix
    #       into a numeric matrix.

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))
    #   dd = drawdowns(x)
    #   par(mfrow = c(3,3)); assetsSeriesPlot(dd); par(mfrow = c(1,1))
    
    # FUNCTION:

    # Settings:
    n = ncol(x)
    if (length(col) == 1) col = rep(col, times = n)

    # Plot:
    seriesPlot(x, ylab = "Series", col = col, ...)

    # Return Value:
    invisible()
}


################################################################################

