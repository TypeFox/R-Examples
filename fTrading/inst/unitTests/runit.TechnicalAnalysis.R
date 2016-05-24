
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
# FUNCTION:                 UTILITY FUNCTIONS:
#  emaTA                     Exponential Moving Average
#  biasTA                    EMA-Bias
#  medpriceTA                Median Price
#  typicalpriceTA            Typical Price
#  wcloseTA                  Weighted Close Price
#  rocTA                     Rate of Change
#  oscTA                     EMA-Oscillator
# FUNCTION:                 OSCILLATOR INDICATORS:
#  momTA                     Momentum
#  macdTA                    MACD
#  cdsTA                     MACD Signal Line
#  cdoTA                     MACD Oscillator
#  vohlTA                    High/Low Volatility
#  vorTA                     Volatility Ratio
# FUNCTION:                 STOCHASTICS INDICATORS:
#  stochasticTA              Stochastics %K/%D, fast/slow
#  fpkTA                     Fast Percent %K
#  fpdTA                     Fast Percent %D
#  spdTA                     Slow Percent %D
#  apdTA                     Averaged Percent %D
#  wprTA                     Williams Percent %R
#  rsiTA                     Relative Strength Index
# FUNCTION:                 DESCRIPTION - MORE INDICATORS:
#  accelTA                   Acceleration
#  adiTA                     AD Indicator
#  adoscillatorTA            AD Oscillator
#  bollingerTA               Bollinger Bands
#  chaikinoTA                Chaikin Oscillator
#  chaikinvTA                Chaikin Volatility
#  garmanklassTA             Garman-Klass Volatility
#  nviTA                     Negative Volume Index
#  obvTA                     On Balance Volume
#  pviTA                     Positive Volume Index
#  pvtrendTA                 Price-Volume Trend
#  williamsadTA              Williams AD
#  williamsrTA               Williams R%
# FUNCTION:                 SPLUS LIKE MOVING AVERAGES:
#  SMA                       Computes Simple Moving Average
#  EWMA                      Computes Exponentially Weighted  Moving Average
# FUNCTION:                 DESCRIPTION:
#  .dailyTA                  Computes an indicator for technical analysis
# FUNCTION:                 DESCRIPTION:
#  .tradeSignals             Computes trade signals from trading positions
#  .tradeLengths             Computes trade length from trading signals
#  .hitRate                  Computes hit rates from returns and positions
# FUNCTION:                 DESCRIPTION:
#  .emaSlider                EMA Slider
################################################################################


test.utilityFunctions =
function()
{
    #  emaTA - Exponential Moving Average
    #  biasTA - EMA-Bias
    #  medpriceTA - Median Price
    #  typicalpriceTA - Typical Price
    #  wcloseTA - Weighted Close Price
    #  rocTA - Rate of Change
    #  oscTA - EMA-Oscillator

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    # Data Records:
    x = close = X[, "Close"]
    high   = X[, "High"]
    low    = X[, "Low"]
    open   = X[, "Open"]
    volume = X[, "Volume"]

    # Exponential Moving Average:
    TA = emaTA(x, lambda = 0.1, startup = 0)
    dim(TA)
    head(TA)

    # EMA-Bias:
    TA = biasTA(x, lag = 5)
    dim(TA)
    head(TA)

    # Median Price:
    TA = medpriceTA(high, low)
    dim(TA)
    head(TA)

    # Typical Price:
    TA = typicalpriceTA(high, low, close)
    dim(TA)
    head(TA)

    # Weighted Close Price:
    TA = wcloseTA(high, low, close)
    dim(TA)
    head(TA)

    # Rate of Change:
    TA = rocTA(x, lag = 5)
    dim(TA)
    head(TA)

    # EMA-Oscillator:
    TA = oscTA(x, lag1 = 25, lag2 = 65)
    dim(TA)
    head(TA)

    # Return Value
    return()
}


################################################################################


test.oscillatorIndicators =
function()
{
    #  momTA - Momentum
    #  macdTA - MACD
    #  cdsTA - MACD Signal Line
    #  cdoTA - MACD Oscillator
    #  vohlTA - High/Low Volatility
    #  vorTA - Volatility Ratio

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    # Data Records:
    x = close = X[, "Close"]
    high   = X[, "High"]
    low    = X[, "Low"]
    open   = X[, "Open"]
    volume = X[, "Volume"]

    # Momentum:
    TA = momTA(x, lag = 5)
    dim(TA)
    head(TA)

    # MACD:
    TA = macdTA(x, lag1 = 12, lag2 = 26)
    dim(TA)
    head(TA)

    # MACD Signal Line:
    TA = cdsTA(x, lag1 = 12, lag2 = 26, lag3 = 9)
    dim(TA)
    head(TA)

    # MACD Oscillator:
    TA = cdoTA(x, lag1 = 12, lag2 = 26, lag3 = 9)
    dim(TA)
    head(TA)

    # High/Low Volatility:
    TA = vohlTA(high, low)
    dim(TA)
    head(TA)

    # Volatility Ratio:
    TA = vorTA(high, low)
    dim(TA)
    head(TA)

    # Return Value:
    return()
}


################################################################################


test.stochasticsIndicators =
function()
{
    # stochasticTA - Stochastics %K/%D, fast/slow
    # fpkTA - Fast Percent %K
    # fpdTA - Fast Percent %D
    # spdTA - Slow Percent %D
    # apdTA - Averaged Percent %D
    # wprTA - Williams Percent %R
    # rsiTA - Relative Strength Index

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    # Data Records:
    x = close = X[, "Close"]
    high   = X[, "High"]
    low    = X[, "Low"]
    open   = X[, "Open"]
    volume = X[, "Volume"]

    # Fast Stochstic:
    # Note, returns a 2-colum series as output ...
    TA = stochasticTA(close, high, low, lag1 = 5, lag2 = 3, type = "fast")
    dim(TA)
    head(TA, 10)

    # Slow Stochstic:
    # Note, returns a 2-colum series as output ...
    TA = stochasticTA(close, high, low, lag1 = 5, lag2 = 3, lag3 = 5,
        type = "slow")
    dim(TA)
    head(TA, 10)

    # Fast Percent K:
    TA = fpkTA(close, high, low, lag = 5)
    dim(TA)
    head(TA,10)

    # Fast Percent D:
    TA = fpdTA(close, high, low, lag1 = 5, lag2 = 3)
    dim(TA)
    head(TA, 10)

    # Slow Percent %D
    TA = spdTA(close, high, low, lag1 = 5, lag2 = 3, lag3 = 9)
    dim(TA)
    head(TA, 10)

    # Averaged Percent %D
    TA = apdTA(close, high, low, lag1 = 5, lag2 = 3, lag3 = 9, lag4 = 9)
    dim(TA)
    head(TA, 10)

    # Williams Percent %R
    TA = wprTA(close, high, low, lag = 5)
    dim(TA)
    head(TA, 10)

    # Relative Strength Index
    TA = rsiTA(close, lag = 14)
    dim(TA)
    head(TA, 10)

    # Return Value:
    return()
}


################################################################################


test.moreIndicators =
function()
{
    # accelTA - Acceleration
    # adiTA - AD Indicator
    # adoscillatorTA - AD Oscillator
    # bollingerTA - Bollinger Bands
    # chaikinoTA - Chaikin Oscillator
    # chaikinvTA - Chaikin Volatility
    # garmanklassTA - Garman-Klass Volatility
    # nviTA - Negative Volume Index
    # obvTA - On Balance Volume
    # pviTA - Positive Volume Index
    # pvtrendTA - Price-Volume Trend
    # williamsadTA - Williams AD
    # williamsrTA- Williams R%

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    x = close = X[, "Close"]
    high   = X[, "High"]
    low    = X[, "Low"]
    open   = X[, "Open"]
    volume = X[, "Volume"]

    # Acceleration
    TA = accelTA(x, n = 3)
    dim(TA)
    head(TA, 10)

    # AD Indicator
    TA = adiTA(high, low, close, volume)
    dim(TA)
    head(TA, 10)

    # AD Oscillator
    TA = adoscillatorTA(open, high, low, close)
    dim(TA)
    head(TA, 10)

    # Bollinger Bands
    TA = bollingerTA(x, lag = 5, n.sd = 2)
    dim(TA)
    head(TA, 10)

    # Chaikin Oscillator
    TA = chaikinoTA(high, low, close, volume, lag1 = 10, lag2 = 3)
    dim(TA)
    head(TA, 10)

    # Chaikin Volatility
    TA = chaikinvTA(high, low, lag1 = 5, lag2 = 5)
    dim(TA)
    head(TA, 10)

    # Garman-Klass Volatility
    TA = garmanklassTA(open, high, low, close)
    dim(TA)
    head(TA, 10)

    # Negative Volume Index
    TA = nviTA(close, volume)
    dim(TA)
    head(TA, 10)

    # On Balance Volume
    TA = obvTA(close, volume)
    dim(TA)
    head(TA, 10)

    # Positive Volume Index
    TA = pviTA(close, volume)
    dim(TA)
    head(TA, 10)

    # Price-Volume Trend
    TA = pvtrendTA(close, volume)
    dim(TA)
    head(TA, 10)

    # Williams AD
    TA = williamsadTA(high, low, close)
    dim(TA)
    head(TA, 10)

    # Williams R%
    TA = williamsrTA(high, low, close, lag = 5)
    dim(TA)
    head(TA, 10)

    # Return Value:
    return()
}


################################################################################


test.splusLikeIndicators =
function()
{
    # SMA - Computes Simple Moving Average
    # EWMA - Computes Exponentially Weighted  Moving Average

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    # Data Records:
    x = close = X[, "Close"]
    high   = X[, "High"]
    low    = X[, "Low"]
    open   = X[, "Open"]
    volume = X[, "Volume"]

    # SMA:
    TA = SMA(x, n = 5)
    dim(TA)
    head(TA)

    # EMA - Using Decay Length:
    TA = EWMA(x, 25)
    dim(TA)
    head(TA)

    # EMA - Using lambda:
    TA = EWMA(x, 2/(25+1))
    dim(TA)
    head(TA)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.dailyTA =
function()
{
    # .dailyTA
    #   Computes an indicator for technical analysis

    # Data from fEcofin:
    X = MSFT
    print(head(X))

    # EMA - Daily TA:
    TA = .dailyTA(X, "ema", select = "Close", lag = 5)
    head(TA)

    # MACD - Daily TA:
    TA = .dailyTA(X, "macd", select = "Close", lag = c(lag1 = 12, lag2 = 26))
    head(TA)

    # ...

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.tradingFunctions =
function()
{
    # .tradeSignals - Computes trade signals from trading positions
    # .tradeLengths - Computes trade length from trading signals
    # .hitRate - Computes hit rates from returns and positions

    # Positions:
    long = +1
    short = -1
    neutral = 0
    tradePositions = c(+1, +1, +1, -1, -1, +1, +1, -1, +1, +1, +1, -1)
    tradeReturns = rnorm(12)

    # Compute Trade Signals:
    Positions = timeSeries(tradePositions, timeCalendar(), units = "Position")
    Positions
    tradeSignals = .tradeSignals(Positions)
    tradeSignals

    # Compute Trade Lengths:
    tradeLengths = .tradeLengths(tradeSignals)
    tradeLengths

    # Compute Hit Rates:
    .hitRate(tradeReturns, tradePositions)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.emaSlider =
function()
{
    .emaSlider =
    function(x)
    {   # A function implemented by Diethelm Wuertz

        # Description
        #   Displays the selected technical indicator

        # FUNCTION:

        # Internal Function:
        refresh.code = function(...)
        {
            # Sliders:
            lambda1  = .sliderMenu(no = 1)
            lambda2  = .sliderMenu(no = 2)
            startup  = .sliderMenu(no = 3)

            # Compute Data:
            seriesPlot(x)
            ema1 = emaTA(x, lambda1, startup)
            N1 = ceiling(2/lambda1)-1
            lines(ema1, col = "red")
            ema2 = emaTA(x, lambda2, startup)
            N2 = ceiling(2/lambda2)-1
            lines(ema2, col = "green")
            mText = paste("EMA1 =", N1, "|", "EMA2 =", N2)
            mtext(mText, side = 4, adj = 0, cex = 0.7, col = "grey")

            # Difference:
            seriesPlot(ema2-ema1, type = "h")
            lines(ema2-ema1, col = "red")

            # Reset Frame:
            par(mfrow = c(2, 1), cex = 0.7)
        }

        # Open Slider Menu:
        N = min(10, dim(x)[1])
        print(N)
        .sliderMenu(refresh.code,
           names =       c(  "lamda1", "lamda2", "startup" ),
           minima =      c(   0.01,     0.01,     0        ),
           maxima =      c(   0.99,     0.99,     N        ),
           resolutions = c(   0.01,     0.01,     1        ),
           starts =      c(   0.10,     0.25,     0        ))
    }

    # Chart:
    # .emaSlider(tS)
    NA

    # Return Value:
    return()
}


################################################################################

