
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port:
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               TAILORED PLOT FUNCTIONS:
#  seriesPlot              Returns a tailored return series plot
#  histPlot                Returns a tailored histogram plot
#  densityPlot             Returns a tailored kernel density estimate plot
#  qqnormPlot              Returns a tailored normal quantile-quantile plot
# FUNCTION:               BASIC STATISTICS:
#  basicStats              Returns a basic statistics summary
# FUNCTION:               DESCRIPTION:
#  .distCheck              Checks consistency of distributions
# FUNCTION:               SPLUS FUNCTIONALITY:
#  stdev                   S-PLUS: Returns the standard deviation of a vector
################################################################################


test.seriesPlot =
function()
{
    # Description:
    #   Returns a tailored return series plot

    # Time Series:
    tS = MSFT[, "Close"]

    # Series Plot:
    par(mfrow = c(2, 1))
    seriesPlot(tS)
    seriesPlot(tS)
    points(tS, col = "red", pch = 19, cex = 0.7)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.histPlot =
function()
{
    # Description:
    #   Returns a tailored histogram plot

    # Time Series:
    tS = MSFT[, "Close"]

    # Histogram Plot:
    par(mfrow = c(1, 1))
    histPlot(tS)
    histPlot(tS, add.fit = FALSE)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.densityPlot =
function()
{
    # Description:
    #   Returns a tailored kernel density estimate plot

    # Time Series:
    tS = MSFT[, "Close"]

    # Density Plot:
    par(mfrow = c(1, 1))
    densityPlot(tS)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.qqnormPlot =
function()
{
    # Description:
    #   Returns a tailored normal quantile-quantile plot

    # Time Series:
    tS = MSFT[, "Close"]

    # Quantile Plot:
    par(mfrow = c(1, 1))
    qqnormPlot(tS)

    # Return Value:
    return()
}



# ------------------------------------------------------------------------------


test.basicStats =
function()
{
    # Description:
    #   Returns a basic statistics summary

    # Time Series:
    tS = MSFT

    Close = tS[, "Close"]

    # Univariate timeSeries - basicStats(x, ci = 0.95)
    basicStats(Close)
    # basicStats(as.numeric(Close))
    basicStats(as.matrix(Close))
    basicStats(as.data.frame(Close))
    basicStats(as.ts(Close))

    # Multivariate - timeSeries - basicStats(x, ci = 0.95)
    basicStats(tS)
    basicStats(as.matrix(tS))
    basicStats(as.data.frame(tS))
    basicStats(as.ts(tS))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.distCheck =
function()
{
    # Description:
    #   Checks consistency of distributions

    # Arguments:
    #   .distCheck(fun = "norm", n = 1000, seed = 4711, ...)

    # Normal Distribution Check:
    fBasics:::.distCheck()

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.stdev =
function()
{
    # Description:
    #   S-PLUS: Returns the standard deviation of a vector

    # Time Series:
    tS = MSFT

    # stdev - Univariate:
    tU = tS[, 1]

    # S-Plus Compatible:
    stdev(tU)
    ## stdev(as.numeric(tU))                                        ## CHECK !!!
    stdev(as.vector(tU))
    stdev(as.ts(tU))

    # Base R:
    sd(tU)
    ## sd(as.numeric(tU))                                           ## CHECK !!!
    sd(as.vector(tU))
    sd(as.ts(tU))


    # Return Value:
    return()
}


################################################################################

