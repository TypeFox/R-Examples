
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
#  assetsStarsPlot             Draws segment/star diagrams of asset sets
# FUNCTION:                   DESCRIPTION:
#  assetsBasicStatsPlot        Displays a segment plot of basic return stats
#  assetsMomentsPlot           Displays a segment plot of distribution moments
#  assetsBoxStatsPlot          Displays a segment plot of box plot statistics
#  assetsNIGFitPlot            Displays a segment plot NIG parameter estimates
################################################################################


assetsStarsPlot <-
    function(x, method = c("segments", "stars"), locOffset = c(0, 0),
    keyOffset = c(0, 0), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Draws segment or star diagrams of a multivariate data set.

    # Arguments
    #   x - a numeric feature matrix of assets. Each column represents
    #       an individual asset.

    # Example:
    #   x = as.timeSeries(data(LPP2005REC))
    #   X = basicStats(x)[-(1:2), 1:6]
    #   assetsStarsPlot(X, main = "Basic Statistics", keyOffset = -0.5)

    # FUNCTION:

    # Settings:
    method = match.arg(method)
    if (method == "segments") draw.segments = TRUE else draw.segments = FALSE

    # Compute Locations:
    xCol = ncol(x)
    yCol = nrow(x)
    NY = NX = ceiling(sqrt(xCol))

    if (NX*NY == xCol) NY = NY + 1

    loc = NULL
    for (nx in 1:NY)
        for (ny in 1:NX)
            loc = rbind(loc, c(nx, ny))
    loc = loc[1:xCol, ]
    loc[, 2] = NY + 1 - loc[, 2]

    loc[, 1] = loc[, 1] - locOffset[1]
    loc[, 2] = loc[, 2] - locOffset[2]

    # Stars:
    palette(rainbow(12, s = 0.6, v = 0.75))
    ans = stars(t(x), mar = c(0,0,0,0), #mar = c(4, 2.8, 2.8, 4),
        locations = loc,
        len = 0.4,
        xlim = c(1, NX+0.5),
        ylim = c(0, NY+1),
        key.loc = c(NX + 1, 1) + keyOffset,
        draw.segments = draw.segments, ... )
    # box()

    # Return Value:
    invisible(ans)
}


# ------------------------------------------------------------------------------


assetsBasicStatsPlot <-
    function(x, par = TRUE, oma = c(0,0,0,0), mar = c(4, 4, 4, 4),
    keyOffset = c(-0.65, -0.50), main = "Assets Statistics",
    title = "Assets", titlePosition = c(3, 3.65),
    description = "Basic Returns Statistics", descriptionPosition = c(3, 3.50),
    ...)
{
    # A function Implemented by Diethelm Wuertz

    # Description:
    #   Displays a segment plot of basic return statistics

    # Note:
    #    The Default Settings are made for a portfolio with
    #       7 to 9 assets.

    # FUNCTION:

    # Plot:
    if (par) par(mfrow = c(1, 1), oma = oma, mar = mar)
    X = basicStats(x)[-(1:2), ]
    assetsStarsPlot(X, keyOffset = keyOffset, ...)
    text(titlePosition[1], titlePosition[2], adj = 0,
        title, cex = 1.25)
    text(descriptionPosition[1], descriptionPosition[2], adj = 0,
        description, cex = 1.1)
    title(main = main)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


assetsMomentsPlot <-
    function(x, par = TRUE, oma = c(0,0,0,0), mar = c(4, 4, 4, 4),
    keyOffset = c(-0.65, -0.50), main = "Assets Statistics",
    title = "Assets", titlePosition = c(3, 3.65),
    description = "Moments Statistics", descriptionPosition = c(3, 3.50), 
    ...)
{
    # A function Implemented by Diethelm Wuertz

    # Description:
    #   Displays a segment plot of distribution moments

    # Note:
    #    The Default Settings are made for a portfolio with
    #       7 to 9 assets.

    # FUNCTION:

    # Plot:
    if(par) par(mfrow = c(1, 1), oma = oma, mar = mar)
    param = NULL
    for (i in 1:dim(x)[2]) {
        X = as.vector(series(x[, i]))
        fit = c(mean = mean(X), stdev = sd(X),
            skewness = skewness(X), kurtosis = kurtosis(X))
        param = cbind(param, fit)
    }
    colnames(param) = colnames(x)
    assetsStarsPlot(param, keyOffset = keyOffset, ...)
    text(titlePosition[1], titlePosition[2], adj = 0,
        title, cex = 1.25)
    text(descriptionPosition[1], descriptionPosition[2], adj = 0,
        description, cex = 1.1)
    title(main = main)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


assetsBoxStatsPlot <-
    function(x, par = TRUE, oma = c(0,0,0,0), mar = c(4, 4, 4, 4),
    keyOffset = c(-0.65, -0.50), main = "Assets Statistics",
    title = "Assets", titlePosition = c(3, 3.65),
    description = "Box Plot Statistics", descriptionPosition = c(3, 3.50), 
    ...)
{
    # A function Implemented by Diethelm Wuertz

    # Description:
    #   Displays a segment plot of box plot statistics

    # Note:
    #    The Default Settings are made for a portfolio with
    #       7 to 9 assets.

    # FUNCTION:

    # Plot:
    if(par) par(mfrow = c(1, 1), oma = oma, mar = mar)
    bp = assetsBoxPlot(x, plot = FALSE)
    ans = assetsStarsPlot(abs(bp$stats), keyOffset = keyOffset, ...)
    text(titlePosition[1], titlePosition[2], adj = 0,
        title, cex = 1.25)
    text(descriptionPosition[1], descriptionPosition[2], adj = 0,
        description, cex = 1.1)
    title(main = main)

    # Return Value:
    invisible(ans)
}


# ------------------------------------------------------------------------------


assetsNIGFitPlot <-
    function(x, par = TRUE, oma = c(0,0,0,0), mar = c(4, 4, 4, 4),
    keyOffset = c(-0.65, -0.50), main = "Assets Statistics",
    title = "Assets", titlePosition = c(3, 3.65),
    description = "NIG  Parameters", descriptionPosition = c(3, 3.50), ...)
{
    # A function Implemented by Diethelm Wuertz

    # Description:
    #   Displays a segment plot NIG parameter estimates

    # Note:
    #    The Default Settings are made for a portfolio with
    #       7 to 9 assets.

    # FUNCTION:

    # Plot:
    param = NULL
    for (i in 1:dim(x)[2]) {
        fit = nigFit(x[, i], doplot = FALSE, trace = FALSE)
        param = cbind(param, fit@fit$estimate)
    }
    if(par) par(mfrow = c(1, 1), oma = oma, mar = mar)
    colnames(param) = colnames(x)
    rownames(param) = c("alpha", "beta", "delta", "mu")
    assetsStarsPlot(param, keyOffset = keyOffset)
    text(titlePosition[1], titlePosition[2], adj = 0,
        title, cex = 1.25)
    text(descriptionPosition[1], descriptionPosition[2], adj = 0,
        description, cex = 1.1)
        title(main = main)

    # Return Value:
    invisible()
}


################################################################################

