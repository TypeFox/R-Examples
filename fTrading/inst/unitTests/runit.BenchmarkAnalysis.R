
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
# FUNCTION:                 BENCHMARK ANALYSIS FUNCTIONS:
#  getReturns                Computes return series given a price series
# FUNCTION:                 DRAWDOWNS:
#  maxDrawDown               Computes the maximum drawdown
# FUNCTION:                 PERFORMANCE RATIOS:
#  sharpeRatio               Calculates the Sharpe Ratio
#  sterlingRatio             Calculates the Sterling Ratio
# FUNCTION:                 OHLC PLOT:
#  ohlcPlot                  Creates a Open-High-Low-Close plot
################################################################################


test.getReturns =
function()
{
    # getReturns - Computes return series given a price series

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    # Get Returns:
    R = getReturns(X)
    head(R)

    # Get Returns:
    R = getReturns(X, percentage = TRUE)
    head(R)

    # Return Value:
    return()
}


################################################################################


test.maxDrawDown =
function()
{
    # maxDrawDown - Computes the maximum drawdown

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    # Closing Prices:
    Close = as.timeSeries(X)[, "Close"]

    # Maximum Draw Down:
    maxDrawDown(Close)

    # Plot:
    plot(Close, type = "l")
    abline(v = as.POSIXct("2000-11-09"), lty = 3, col = "red")
    abline(v = as.POSIXct("2000-12-20"), lty = 3, col = "red")

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.sharpeRatio =
function()
{
    # sharpeRatio - Calculates the Sharpe Ratio

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    # Get Returns:
    R = getReturns(X)

    # Sharpe Ratio:
    sharpeRatio(R[, "Close"])

    # Return Value:
    return()
}


################################################################################


test.sterlingRatio =
function()
{
    # sterlingRatio - Calculates the Sterling Ratio

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    # Get Returns:
    R = getReturns(X)

    # Sterling Ratio:
    sterlingRatio(R[, "Close"])

    # Return Value:
    return()
}


################################################################################


test.ohlcPlot =
function()
{
    #  ohlcPlot - Creates a Open-High-Low-Close plot

    # Data from fEcofin:
    myFinCenter <<- "GMT"
    X = MSFT
    print(head(X))

    # Get Returns:
    R = returns(X)[, -5]
    Y = alignDailySeries(X, method = "fillNA", include.weekends = TRUE)

    # Plot:
    # ohlcPlot(as.ts(R))                                             # CHECK !!!

    # Return Value:
    return()
}


################################################################################

