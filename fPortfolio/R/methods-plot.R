
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA


################################################################################
# FUNCTION:                     DESCRIPTION:
#  plot.fPORTFOLIO               S3 Plot method for 'fPORTFOLIO' objects
# FUNCTION:                     DESCRIPTION:
#  .fPortfolio.plot1..8          Internal plot functions
################################################################################


plot.fPORTFOLIO <-
    function(x, which = "ask", control = list(), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class 'fPORTFOLIO'

    # FUNCTION:

    # Control Parameters:
    Statistics = getStatistics(x)

    # Use default, if xlim and ylim is not specified ...
    mu = Statistics$mu
    Sigma = Statistics$Sigma
    N = length(mu)
    yLim = range(mu) + 0.25*c(-diff(range(mu)), diff(range(mu)))

    # First, take care that all assets appear on the plot ...
    # sqrtSig = sqrt(diag(Sigma))
    # xLimAssets = c(
    #    min(sqrtSig),
    #    max(sqrtSig))+ c(-0.4*diff(range(sqrtSig)), 0.1*diff(range(sqrtSig)))
    xRange = range(frontierPoints(x)[, 1])
    xDiff = diff(xRange)
    xLimAssets = c(xRange[1] - 2.5*xDiff/10, xRange[2] + xDiff/10)

    # ... second take care that the whole frontier appears on the plot:
    fullFrontier = frontierPoints(x)
    xLimFrontier = range(fullFrontier[, 1])
    xLim = range(c(xLimAssets, xLimFrontier))

    # Control List:
    # YC: merge with user control list
    con <- c(control, frontierPlotControl())
    # YC: remove double entries and keep user args
    con <- con[unique(names(con))]
    attr(x, "control") <- con

    # Plot:
    interactivePlot(
        x,
        choices = c(
            "Plot Efficient Frontier",
            "Add Minimum Risk Portfolio",
            "Add Tangency Portfolio",
            "Add Risk/Return of Single Assets",
            "Add Equal Weights Portfolio",
            "Add Two Asset Frontiers [LongOnly Only]",
            "Add Monte Carlo Portfolios",
            "Add Sharpe Ratio [Markowitz PF Only]"),
        plotFUN = c(
            ".fportfolio.plot.1", ".fportfolio.plot.2", ".fportfolio.plot.3",
            ".fportfolio.plot.4", ".fportfolio.plot.5", ".fportfolio.plot.6", 
            ".fportfolio.plot.7", ".fportfolio.plot.8"),
        which = which)

    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.fportfolio.plot.1 <-
    function(x)
{
    # Description:
    #   Plot Efficient Frontier

    # FUNCTION:

    # Control:
    con = attr(x, "control")
    Type = getType(x)
    if (Type == "MV") {
        xLab = "Mean-Var Target Risk"
    } else if (Type == "CVaR") {
        xLab = "-CVaR Target Risk"
    }

    # Plot:
    frontierPlot(
        object = x, xlim = con$xlim, ylim = con$ylim, 
        pch = 19, cex = 0.75, title = FALSE,
        las = ifelse(is.null(con$las), 0, con$las))
    title(
        main = ifelse(is.null(con$main), "Efficient Frontier", con$main),
        xlab = ifelse(is.null(con$xlab), xLab, con$xlab),
        ylab = ifelse(is.null(con$ylab), "Target Return", con$ylab))

}


# ------------------------------------------------------------------------------


.fportfolio.plot.2 <-
    function(x)
{
    # Description:
    #   Add Minimum Risk Portfolio

    # FUNCTION:

    # Control:
    con = attr(x, "control")

    # Plot:
    minvariancePoints(
        object = x,
        col = con$minvariance.col,
        cex = con$minvariance.cex,
        pch = 19)
}


# ------------------------------------------------------------------------------


.fportfolio.plot.3 <-
    function(x)
{
    # Description:
    #   Add Tangency Portfolio

    # FUNCTION:

    # Control:
    con = attr(x, "control")

    # Plot:
    tangencyPoints(
        object = x,
        col = con$tangency.col,
        cex = con$tangency.cex,
        pch = 17)
    tangencyLines(
        object = x,
        col = con$tangency.col,
        cex = con$tangency.cex)
}


# ------------------------------------------------------------------------------


.fportfolio.plot.4 <-
    function(x)
{
    # Description:
    #   Add Risk/Return of Single Assets

    # FUNCTION:

    # Control:
    con = attr(x, "control")

    # Plot:
    Palette = match.fun(con$singleAsset.col)
    col = Palette(getNAssets(x))
    singleAssetPoints(
        object = x,
        col = col,
        cex = con$singleAsset.cex,
        pch = 18)
}


# ------------------------------------------------------------------------------


.fportfolio.plot.5 <-
    function(x)
{
    # Description:
    #   Add Equal Weights Portfolio

    # FUNCTION:

    # Control:
    con = attr(x, "control")

    # Plot:
    equalWeightsPoints(
        object = x,
        col = con$equalWeights.col,
        cex = con$equalWeights.cex,
        pch = 15)
}


# ------------------------------------------------------------------------------


.fportfolio.plot.6 <-
    function(x)
{
    # Description:
    #   Add Two Asset Frontiers [0-1 PF Only]

    # FUNCTION:

    # Control:
    con = attr(x, "control")

    # Lines:
    lines(frontierPoints(object = x), col = "grey")
    twoAssetsLines(object = x, col = con$twoAssets.col)

    # Points:
    Palette = match.fun(con$singleAsset.col)
    col = Palette(getNAssets(x))
    singleAssetPoints(
        object = x,
        col = col,
        cex = con$singleAsset.cex,
        pch = 18)
}


# ------------------------------------------------------------------------------


.fportfolio.plot.7 <-
    function(x)
{
    # Description:
    #   Add Monte Carlo Portfolios

    # FUNCTION:

    # Control:
    con = attr(x, "control")

    # Plot:
    monteCarloPoints(
        object = x,
        col = con$monteCarlo.col,
        cex = con$monteCarlo.cex,
        mcSteps = con$mcSteps)
}


# ------------------------------------------------------------------------------


.fportfolio.plot.8 <-
    function(x)
{
    # Description:
    #   Add Sharpe Ratio [MV PF Only]
    # FUNCTION:

    # Control:
    con = attr(x, "control")

    # Plot:
    sharpeRatioLines(
        object = x,
        col = con$sharpeRatio.col,
        cex = con$sharpeRatio.cex,
        lty = 3)
}


################################################################################

