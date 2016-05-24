
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
# You should have received a copy of the GNU Library General
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
# FUNCTION          DATA PREPROCESSING:
#  blockMaxima       Returns block maxima from a time series
#  findThreshold     Upper threshold for a given number of extremes
#  pointProcess      Returns peaks over a threshold from a time series
#  deCluster         Declusters a point process
################################################################################


test.blockMaxima =
function()
{
    # blockMaxima - Returns block maxima from a time series

    # blockMaxima(x, block = c("monthly", "quarterly"), doplot = FALSE)

    # Time Series Data:
    x = MSFT[, "Close"]
    x.ret = 100*returns(x)
    head(x.ret)
    class(x.ret)

    # Monthly Block Maxima:
    ans = blockMaxima(x.ret, block = "monthly", doplot = TRUE)
    print(ans)

    # Quarterly Block Maxima:
    ans = blockMaxima(x.ret, block = "q", doplot = TRUE)
    print(ans)

    # 20-Days Block Maxima:
    ans = blockMaxima(x.ret, block = 20, doplot = TRUE)
    print(ans)

    # Numerical Data Vector:
    x.ret = as.vector(x.ret)
    head(x.ret)
    ans = blockMaxima(x.ret, block = 20, doplot = TRUE)
    print(ans)

    # Stops by stopifnot() - Check:
    # blockMaxima(x.ret, block = "month", doplot = TRUE)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.findThreshold =
function()
{
    # findThreshold - Upper threshold for a given number of extremes

    # findThreshold(x, n = floor(0.05*length(as.vector(x))), doplot = FALSE)

    # Time Series Data:
    x = MSFT[, "Close"]
    x.ret = 100*returns(x)
    head(x.ret)
    class(x.ret)

    # Find 99% Threshold:
    par(mfrow = c(2, 2), cex = 0.7)
    par(ask = FALSE)
    findThreshold(x.ret, n = floor(0.01*length(as.vector(x))), doplot = TRUE)

    # Remark - Alternative use ...
    quantile(x.ret, probs = 1 - 0.01)
    quantile(x.ret, probs = 1 - 0.01, type = 1)

    # Find 95% Threshold:
    findThreshold(x.ret, doplot = TRUE)

    # Find 90% Threshold:
    findThreshold(x.ret, n = floor(0.1*length(as.vector(x))), doplot = TRUE)

    # Try if x is a numeric vector:
    findThreshold(as.vector(x.ret), doplot = TRUE)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.pointProcess =
function()
{
    # pointProcess - Returns peaks over a threshold from a time series

    # pointProcess(x, u = quantile(x, 0.95), doplot = FALSE)

    # Time Series Data:
    x = MSFT[, "Close"]
    x.ret = 100*returns(x)
    head(x.ret)
    class(x.ret)

    # Plot Series:
    par(mfrow = c(2, 1), cex = 0.7)
    par(ask = FALSE)
    # plot(x.ret, type = "l", main = "Series")
    # abline(h = 0, col = "red", lty = 3)
    # or use ...
    seriesPlot(x.ret)

    # Point Process:
    pp = pointProcess(x.ret, u = quantile(x.ret, 0.8))
    pp
    plot(pp, type = "b", main = "Point Process")
    abline(h = 0, col = "red", lty = 3)
    # Try seriesPlot(pp)
    # ... add points in graph

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.deCluster =
function()
{
    # deCluster - Declusters a point process

    # deCluster(x, run = 20, doplot = TRUE)

    # Time Series Data:
    x = MSFT[, "Close"]
    x.ret = 100*returns(x)
    head(x.ret)
    class(x.ret)

    # Decluster Time Series:
    tS = deCluster(x = x.ret, run = 3)
    print(tS)

    # Return Value:
    return()
}


################################################################################

