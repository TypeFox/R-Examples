
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
#  assetsQQNormPlot            Displays normal qq-plots of individual assets
################################################################################


assetsQQNormPlot =
    function(x, col = "steelblue", skipZeros = FALSE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays normal qq-plots of individual assets

    # Arguments:
    #   x - a timeSeries object or any other rectangular object
    #       which can be transformed by the function as. matrix
    #       into a numeric matrix.
    #   which - an integer value or vector specifying the number(s)
    #       of the assets which are selected to be plotted.

    # FUNCTION:

    # Settings:
    n = ncol(x)
    if (length(col) == 1) col = rep(col, times = n)

    # Plot:
    for (i in 1:n) {
        X = x[, i]
        if (skipZeros) X = X[series(X) != 0]
        qqnormPlot(X, col = col[i], ...)
    }

    # Return Value:
    invisible()
}


################################################################################


assetsHistPairsPlot <-
    function(x, bins = 30, method = c("square", "hex"), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays bivariate Histogram Plot

    # FUNCTION:

    # Match Arguments:
    method = match.arg(method)

    # Check:
    stopifnot(ncol(x) == 2)

    # Histogram Plot:
    X = as.vector(x[, 1])
    Y = as.vector(x[, 2])
    if (method == "square") {
        ans = squareBinning(x = X, y= Y, bins = bins)
    } else if (method == "hex") {
        ans = hexBinning(x = X, y = Y, bins = bins)
    }

    # Plot:
    plot(ans, ...)

    # Return Value:
    invisible(ans)
}


################################################################################

