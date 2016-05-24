
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
#  assetsHistPairsPlot         Displays a bivariate histogram plot 
################################################################################


assetsHistPairsPlot <-
    function(x, bins = 30, method = c("square", "hex"), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays bivariate Histogram Plot
    
    # Arguments:
    #   x - timeSeries
    #   bins - histogram bins
    #   method - plot method
     
    # Example:
    #   x <- 100 * as.timeSeries(data(LPP2005REC))[, c("SBI", "SPI")]
    #   assetsHistPairsPlot(x, bins = 20)
    #   assetsHistPairsPlot(x, bins = 20, method = "hex")
    
    # FUNCTION:

    # Match Arguments:
    method <- match.arg(method)

    # Check:
    stopifnot(ncol(x) == 2)

    # Histogram Plot:
    X <- as.vector(x[, 1])
    Y <- as.vector(x[, 2])
    if (method == "square") {
        ans <- squareBinning(x = X, y= Y, bins = bins)
    } else if (method == "hex") {
        ans <- hexBinning(x = X, y = Y, bins = bins)
    }

    # Plot:
    plot(ans, ...)

    # Return Value:
    invisible(ans)
}


################################################################################

