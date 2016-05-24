
# Rmetrics is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# Rmetrics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################


test.timeSeries =
function()
{
    #  timeSeries - Creates a 'timeSeries' object from scratch

    # Settings:
    setRmetricsOptions(myFinCenter = "GMT")
    set.seed(4711)
    data = matrix(round(rnorm(12), 3))
    data
    class(data)
    charvec = format(timeCalendar(2006))
    charvec
    class(charvec)

    # Compose Univariate daily random sequence
    setRmetricsOptions(myFinCenter = "GMT")
    uTS = timeSeries(data, charvec, units = "uTS")
    series(uTS)
    print(uTS)

    # FinCenter Functionality:
    timeSeries(data, charvec, units = "uTS", zone = "GMT", FinCenter = "GMT")
    timeSeries(data, charvec, units = "uTS", zone = "Zurich", FinCenter = "Zurich")
    timeSeries(data, charvec, units = "uTS", zone = "GMT", FinCenter = "Zurich")
    timeSeries(data, charvec, units = "uTS", zone = "Zurich", FinCenter = "GMT")

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.readSeries =
function()
{
    #  readSeries - Reads from a spreadsheet and creates a 'timeSeries'

    # Load Microsoft Data:
    data(MSFT)
    MSFT.df = as.data.frame(MSFT)

    # Read Data Frame:
    write.table(MSFT.df, file = "msft.dat.csv", sep = ";")
    read.table("msft.dat.csv", sep = ";")

    # Read Time Series:
    # X = readSeries("msft.dat.csv")
    # X = X[1:12, ]
    # class(X)

    # Show Part of Series:
    # head(X)[, 1:5]
    # head(X[, 1:5])
    # head(X[, 1:5], 2)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.returns =
function()
{
    #  returns - Computes returns from a 'timeSeries' object

    # Load Time Series:
    X = MSFT
    head(X)

    # returns :
    OPEN = X[, 1]
    print(OPEN)
    MSFT.RET = returns(OPEN)
    print(MSFT.RET)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.applySeries =
function()
{
    #  applySeries - Applies a function to blocks of a 'timeSeries'

    NA

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.orderStatistics =
function()
{
    #  orderStatistics - Compute order statistic of a 'timeSeries'

    # Load Data:
    X = MSFT
    head(X)

    # returns:
    OPEN = X[, 1]
    print(OPEN)

    # ORDER STATISTICS:
    orderStatistics(OPEN)
    orderStatistics(X[, -5])
    orderStatistics(X[, -5])$Open

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.series =
function()
{
    #  series - Extracts data slot from 'timeSeries' object

    # Load Microsoft Data:
    X = MSFT
    X = X[1:12, ]
    class(X)

    # Return Series:
    OPEN = X[, 1]
    OPEN
    returns(OPEN)

    # Volatility Series:
    abs(returns(OPEN))

    # Data Matrix:
    series(OPEN)
    Y = series(X)
    Y
    class(Y)

    # Position Vector:
    PO = time(OPEN)
    PO
    PX = time(X)
    PX
    class(PX)
    checkEquals(
        target = sum(as.integer(PO - PX)),
        current = 0)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.isUnivariate =
function()
{
    #  isUnivariate     Tests if an object of class 'timeSeries' is univariate

    # Load Microsoft Data:
    X = MSFT
    OPEN = X[, 1]

    # Is Univariate?
    checkTrue(!isUnivariate(X))
    checkTrue(isUnivariate(OPEN))


    checkTrue(isMultivariate(X))
    checkTrue(!isMultivariate(OPEN))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.isMultivariate =
function()
{
    #  isMultivariate - Tests if an object of class 'timeSeries' is multivariate

    # Load Microsoft Data:
    X = MSFT
    OPEN = X[, 1]

    # Is Multivariate?
    checkTrue(isMultivariate(X))
    checkTrue(!isMultivariate(OPEN))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.displayMethods =
function()
{
    #  print.timeSeries     Print method for a 'timeSeries' object
    #  plot.timeSeries      Plot method for a 'timeSeries' object
    #  lines.timeSeries     Lines method for a 'timeSeries' object
    #  points.timeSeries    Points method for a 'timeSeries' object

## FIXME(MM) - if we store this -- make it a package data set!
##     Microsoft Data:
##     MSFT.df = data.frame(matrix(c(
##     20010326, 57.1250, 57.5000, 55.5625, 56.0625,  31559300,
##     20010327, 56.0625, 58.5625, 55.8750, 58.2500,  47567800,
##     20010328, 57.3750, 57.9375, 55.3750, 55.5625,  39340800,
##     20010329, 55.3750, 57.1875, 54.5625, 55.3750,  43492500,
##     20010330, 55.7500, 56.1875, 53.8750, 54.6875,  45600800,
##     20010402, 54.8125, 56.9375, 54.6250, 55.8125,  37962000,
##     20010403, 55.3125, 55.3125, 52.7500, 53.3750,  47093800,
##     20010404, 53.3750, 55.0000, 51.0625, 51.9375,  52023300,
##     20010405, 53.7500, 57.3750, 53.5000, 56.7500,  56682000,
##     20010406, 56.3750, 57.1875, 55.0625, 56.1875,  46311000,
##     20010409, 56.5700, 57.4200, 55.6600, 57.1500,  28147800,
##     20010410, 57.9500, 60.0900, 57.7800, 59.6800,  54599700,
##     20010411, 60.6500, 61.5000, 59.7000, 60.0400,  54939800,
##     20010412, 59.5600, 62.3100, 59.3500, 62.1800,  43760000,
##     20010416, 61.4000, 61.5800, 60.1200, 60.7900,  32928700,
##     20010417, 60.5200, 62.1100, 60.0400, 61.4800,  42574600,
##     20010418, 63.3900, 66.3100, 63.0000, 65.4300,  78348200,
##     20010419, 65.8100, 69.0000, 65.7500, 68.0400,  79687800,
##     20010420, 70.3000, 71.1000, 68.5000, 69.0000,  96459800,
##     20010423, 68.1100, 68.4700, 66.9000, 68.2500,  46085600,
##     20010424, 68.2000, 69.9300, 67.1400, 67.5500,  44588300,
##     20010425, 67.5700, 69.7900, 67.2500, 69.6900,  38372000,
##     20010426, 70.0700, 71.0000, 68.2500, 69.1300,  59368800,
##     20010427, 69.5300, 69.6800, 66.2100, 67.1200,  60786200,
##     20010430, 68.5300, 69.0600, 67.6800, 67.7500,  37184100,
##     20010501, 67.6600, 70.3000, 67.6000, 70.1700,  41851400,
##     20010502, 71.0000, 71.1500, 69.3500, 69.7600,  46432200,
##     20010503, 69.2500, 70.1800, 68.1400, 68.5300,  33136700,
##     20010504, 68.0000, 71.0500, 67.9600, 70.7500,  59769200,
##     20010507, 70.8300, 72.1500, 70.7000, 71.3800,  54678100),
##     byrow = TRUE, ncol = 6))
##     colnames(MSFT.df) = c("YYMMDD", "Open", "High", "Low", "Close", "Volume")

    # Data:
    X = MSFT
    X = X[1:12, ]
    OPEN = X[, 1]

    # Print:
    print(X)
    print(OPEN)

    # Plot:
    par(mfrow = c(1, 1))
    plot(OPEN, type = "l")

    # GMT - Plot:
    tC = timeCalendar(2006, 1, 1, 0:23, 0, 0, zone = "GMT", FinCenter = "GMT")
    tS = timeSeries(data = matrix(rnorm(24), ncol = 1), charvec = tC)
    plot(tS)

    # Zurich - Plot:
    tC = timeCalendar(2006, 1, 1, 0:23, 0, 0, zone = "GMT", FinCenter = "Zurich")
    tS = timeSeries(data = matrix(rnorm(24), ncol = 1), charvec = tC,
        zone = "GMT", FinCenter = "Zurich")
    plot(tS)

    # New York - Plot:
    tC = timeCalendar(2006, 1, 1, 0:23, 0, 0, zone = "GMT", FinCenter = "NewYork")
    tS = timeSeries(data = matrix(rnorm(24), ncol = 1), charvec = tC,
        zone = "GMT", FinCenter = "NewYork")
    plot(tS, type = "h")
    lines (tS, col = "red",  lty = 3)
    points(tS, col = "blue", pch = 19)
    abline(h=0, col = "grey")

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.dummyDailySeries =
function()
{
    #  dummyDailySeries - Creates a dummy daily 'timeSeries' object

    # Create Dummy Time Series:
    setRmetricsOptions(myFinCenter = "GMT")
    tS = dummyDailySeries(matrix(rnorm(12)))
    print(tS)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.alignDailySeries =
function()
{
    # alignDailySeries - Aligns a 'timeSeries' object to new positions

    # Time Series:
    setRmetricsOptions(myFinCenter = "GMT")
    tS = MSFT[1:25, ]
    print(tS)
    dim(tS)

    # Align Daily Series:
    alignDailySeries(tS, method = "interp")

    # Align Daily Series:
    alignDailySeries(tS, method = "fillNA")

    # Align Daily Series:
    alignDailySeries(tS, method = "fillNA", include.weekends = TRUE)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


## DW >
## test.ohlcDailyPlot =
## function()
## {
##     # ohlcDailyPlot - Plots open–high–low–close bar chart
## 
##     # Price or Incdex Series:
##     setRmetricsOptions(myFinCenter = "GMT")
##     tS = MSFT[1:25, ]
##     print(tS)
##     dim(tS)
##     colnames(tS)
## 
##     # Graph Frame:
##     par(mfrow = c(2, 1), cex = 0.7)
##     ohlcDailyPlot(tS)
## 
##     # Return Value:
##     return()
## }


# ------------------------------------------------------------------------------


test.modelSeries =
function()
{
    if (FALSE) {

        # Move to fArma ...

        # Undocumented Material:
        Matrix = cbind(X = rnorm(10), Y = rnorm(10))
        Matrix = cbind(Matrix, Z = Matrix[, "Y"] - Matrix[, "X"])
        TS = dummyDailySeries(Matrix, units = c("X", "Y", "Z") )
        head(TS)

        .modelSeries(Y ~ ar(2), data = TS, lhs = TRUE)
        .modelSeries(log(abs(Z)) ~ lm(X + sin(Y)), data = TS, fake = TRUE)
        .modelSeries(log(abs(Z)) ~ lm(X + sin(Y)), data = TS, lhs = TRUE)

        .modelSeries(Y ~ ar(2), data = as.data.frame(TS), lhs = TRUE)
        .modelSeries(log(abs(Z)) ~ lm(X + sin(Y)), data = TS, fake = TRUE)
        .modelSeries(log(abs(Z)) ~ lm(X + sin(Y)), data = TS, lhs = TRUE)

        require(timeSeries)
        .modelSeries(Y ~ ar(2), data = rnorm(10))
        .modelSeries(Y ~ ar(2), data = as.ts(rnorm(10)))
        .modelSeries(x ~ arima(2, 0, 1), data = armaSim(n=10))

        .modelSeries(~ ar(2), rnorm(10))

        # attach(TS)                                    # CHECK
        # .modelSeries(Y ~ ar(2), lhs = TRUE)

        .modelSeries(Y ~ ar(2) + garch(1,1), data = rnorm(10))
        .modelSeries(Y ~ ar(2) + garch(1,1), data = rnorm(10), lhs = TRUE)
        .modelSeries(Y ~ ar(2) + garch(1,1), data = TS, lhs = TRUE)

    } else {

        NA

    }

    # Return Value:
    return()
}


################################################################################

