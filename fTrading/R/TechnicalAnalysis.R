
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
################################################################################


# Notations / Data - Prices, Volume, OpenInterest:
#     O    Open         H    High
#     L    Low          C    Close
#     V    Volume       X    one of O H, L, C, V


# ------------------------------------------------------------------------------


emaTA =
function(x, lambda = 0.1, startup = 0)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the Exponential Moving Average Indicator

    # Details:
    #   EXPONENTIAL MOVING AVERAGE:
    #   EMA: EMA(n) = lambda * X(n) + (1-lambda) * EMA(n-1)
    #       lambda = 2 / ( n+1 )

    # Example:
    #   head(emaTA(MSFT[, "Close"]))

    # FUNCTION:

    # Preprocessing:
    TS = is.timeSeries(x)
    y = as.vector(x)

    # EMA:
    if (lambda >= 1) lambda = 2/(lambda+1)
    if (startup == 0) startup = floor(2/lambda)
    if (lambda == 0){
    ema = rep (mean(x),length(x))}
    if (lambda > 0){
        ylam = y * lambda
        ylam[1] = mean(y[1:startup])
    ema = filter(ylam, filter = (1-lambda), method = "rec")}
    ema =  as.vector(ema)

    # Convert to timeSeries object:
    if (TS) {
        ema = matrix(ema)
        colnames(ema) = "EMA"
        rownames(ema) = rownames(x)
        series(x) = ema
    } else {
        x = ema
    }

    # Return Value:
    x
}


# ------------------------------------------------------------------------------


biasTA =
function(x, lag = 5)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the Bias Indiacator

    # Example:
    #   head(biasTA(MSFT[, "Close"]))
    #   head(biasTA(rnorm(30)))

    # Details:
    #   BIAS: (X - EMA) / EMA

    # FUNCTION:

    # BIAS:
    xema = emaTA(x, lag)
    bias = (x - xema)/xema
    if (is.timeSeries(bias)) colnames(bias)<-"BIAS"

    # Return Value:
    bias
}


# ------------------------------------------------------------------------------


medpriceTA  =
function(high, low)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the Middle Price Indicator

    # Example:
    #   head(medpriceTA(MSFT[, "High"], MSFT[, "Low"]))
    #   head(medpriceTA(rnorm(30), rnorm(30)))

    # FUNCTION:

    # MEDPRICE:
    medprice = (high + low) / 2
    if (is.timeSeries(medprice)) colnames(medprice)<-"MEDPRICE"

    # Return Value:
    medprice
}


# ------------------------------------------------------------------------------


typicalpriceTA =
function(high, low, close)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the Typical Price Indicator

    # Example:
    #   head(typicalpriceTA(MSFT[, "High"], MSFT[, "Low"], MSFT[, "Close"]))
    #   head(typicalpriceTA(rnorm(30), rnorm(30), rnorm(30)))

    # FUNCTION:

    # Typical Price
    typicalprice = (high + low + close) / 3
    if (is.timeSeries(typicalprice)) colnames(typicalprice)<-"TYPICALPRICE"

    # Return Value:
    typicalprice
}


# ------------------------------------------------------------------------------


wcloseTA =
function(high, low, close)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns Weighted Close Indicator

    # Example:
    #   head(wcloseTA(MSFT[, "High"], MSFT[, "Low"], MSFT[, "Close"]))
    #   head(wcloseTA(rnorm(30), rnorm(30), rnorm(30)))

    # FUNCTION:

    # Weighted Close:
    wclose = (high + low + 2 * close) / 4
    if (is.timeSeries(wclose)) colnames(wclose)<-"WCLOSE"

    # Return Value:
    wclose
}


# ------------------------------------------------------------------------------


rocTA =
function(x, lag = 5)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns rate of Change Indicator

    # Examples:
    #   head(rocTA(MSFT[, "Close"]))
    #   head(rocTA(rnorm(30)))

    # Details:
    #   RATE OF CHANGE INDICATOR:
    #   ROC: (X(n) - X(n-k) ) / X(n)

    # FUNCTION:

    # Rate of Change:
    if (is.timeSeries(x)) {
        roc = diff(x, lag = lag, pad = 0) / x
        colnames(roc)<-"ROC"
    } else {
        roc = diff(x, lag = lag)
        roc = c(rep(0, times = lag), roc) / x
    }

    # Return Value:
    roc
}


# ******************************************************************************


oscTA =
function(x, lag1 = 25, lag2 = 65)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns EMA Oscillator Indicator

    # Examples:
    #   head(oscTA(MSFT[, "Close"]))
    #   head(oscTA(rnorm(30)))

    # Details:
    #   EMA OSCILLATOR INDICATOR:
    #   OSC: (EMA_LONG - EMA_SHORT) / EMA_SHORT

    # FUNCTION:

    # Oscillator:
    xema1 = emaTA(x, lag1)
    xema2 = emaTA(x, lag2)
    osc = (xema1 - xema2) / xema2
    if (is.timeSeries(osc)) colnames(osc)<-"OSC"

    # Return Value:
    osc
}


# ------------------------------------------------------------------------------


momTA =
function(x, lag = 25)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns Momentum Indicator

    # Examples:
    #   head(momTA(MSFT[, "Close"]))
    #   head(momTA(rnorm(30)))

    # Details:
    #   MOMENTUM INDICATOR:
    #   MOM: X(n) - X(n-lag)

    # FUNCTION:

    # Momentum:
    if (is.timeSeries(x)) {
        mom = diff(x, lag = lag, pad = 0)
        colnames(mom)<-"MOM"
    } else {
        mom = diff(x, lag = lag)
        mom = c(rep(0, times = lag), mom)
    }

    # Return Value:
    mom
}


# ------------------------------------------------------------------------------


macdTA =
function(x, lag1 = 12, lag2 = 26)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns MA Convergence-Divergence Indicator

    # Details
    #   MACD MA CONVERGENCE DIVERGENCE INDICATOR
    #   Fast MACD e.g. lag1=12, lag=26
    #   MCD: (EMA_SHORT - EMA_LONG)

    # FUNCTION:

    # MACD:
    macd = emaTA(x, lag1) - emaTA(x, lag2)
    if (is.timeSeries(x)) colnames(macd)<-"MACD"

    # Return Result:
    macd
}


# ------------------------------------------------------------------------------


cdsTA =
function(x, lag1 = 12, lag2 = 26, lag3 = 9)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns MACD Slow Signal Line Indicator

    # Details:
    #   MACD SLOW SIGNAL LINE INDICATOR: e.g. lag3=9
    #   SIG: EMA(MCD)

    # FUNCTION:

    # CDS:
    cds = emaTA(macdTA(x, lag1, lag2), lag3)
    if (is.timeSeries(x)) colnames(cds)<-"CDS"

    # Return Result:
    cds
}


# ------------------------------------------------------------------------------


cdoTA <-
function(x, lag1 = 12, lag2 = 26, lag3 = 9)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns MA Convergence-Divergence Oscillator Indicator

    # Details:
    #   MACD - MA CONVERGENCE DIVERGENCE OSCILLATOR:
    #   CDO: MACD - SIG

    # FUNCTION:

    # CDO:
    cdo = macdTA(x, lag1 = lag1, lag2 = lag2) - cdsTA(x, lag3 = lag3)
    if(is.timeSeries(x)) colnames(cdo)<-"CDO"

    # Return Value:
    cdo
}


# ------------------------------------------------------------------------------


vohlTA =
function(high, low)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns High Low Volatility Indicator

    # Details:
    #   HIGH LOW VOLATILITY:
    #   VOHL: high - low

    # FUNCTION:

    # VOHL:
    vohl = high - low
    if(is.timeSeries(vohl)) colnames(vohl)<-"VOHL"

    # Return Value:
    vohl
}


# ------------------------------------------------------------------------------


vorTA =
function(high, low)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns Volatility Ratio Indicator

    # Details:
    #   VOLATILITY RATIO:
    #   VOR: (high-low)/low

    # FUNCTION:

    # VOR:
    vor = (high - low) / low
    if(is.timeSeries(vor)) colnames(vor)<-"VOR"

    # Return Value:
    vor
}


# ------------------------------------------------------------------------------


stochasticTA =
function (close, high, low, lag1 = 5, lag2 = 3, lag3 = 5,
type = c("fast", "slow"))
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns Stochastic Indicators

    # Example:
    #   stochasticTA(high, low, close, lag1 = 5, lag2 = 3, "fast")
    #   stochasticTA(high, low, close, lag1 = 5, lag2 = 3, "slow")

    # FUNCTION:

    # Settings:
    TS = is.timeSeries(close)
    if (TS) {
        stopifnot(isUnivariate(close))
        stopifnot(isUnivariate(high))
    }
    type = match.arg(type)

    # Fast:
    K = fpkTA(close, high, low, lag = lag1)
    D = fpdTA(close, high, low, lag1 = lag1, lag2 = lag2)

    # Slow:
    if (type == "slow") {
        K = emaTA(K, lag3)
        D = emaTA(D, lag3)
    }

    # Indicator:
    if (TS) {
        stochastic = cbind(K, D)
        units = c(paste(type, "K", sep = ""), paste(type, "D", sep = ""))
        colnames(stochastic)<-units
    } else {
        stochastic = cbind(K, D)
    }

    # Return Value:
    stochastic
}


# ------------------------------------------------------------------------------


fpkTA =
function(close, high, low, lag = 5)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns Fast %K Indicator

    # FUNCTION:

    # Settings:
    TS = is.timeSeries(close)
    if (TS) {
        X = close
        close = as.vector(close)
        high = as.vector(high)
        low = as.vector(low)
    }

    # Minimum:
    minlag = function(x, lag) {
        xm = x
        for (i in 1:lag) {
            x1 = c(x[1], x[1:(length(x)-1)])
            xm = pmin(xm,x1)
            x = x1
        }
        xm
    }

    # Maximum:
    maxlag = function(x, lag) {
        xm = x
        for (i in 1:lag) {
            x1 = c(x[1], x[1:(length(x)-1)])
            xm = pmax(xm,x1)
            x = x1
       }
       xm
    }

    # Result:
    xmin = minlag(low, lag)
    xmax = maxlag(high, lag)
    fpk = (close - xmin ) / (xmax -xmin)

    # to timeSeries:
    if (TS) {
        fpk = matrix(fpk)
        rownames(fpk) = rownames(series(X))
        colnames(fpk) = "FPK"
        series(X) = fpk
    } else {
        X = fpk
    }

    # Return Value:
    X
}


# ------------------------------------------------------------------------------


fpdTA =
function(close, high, low, lag1 = 5, lag2 = 3)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns Fast %D Indicator

    # Details:
    #   FAST %D INDICATOR: EMA OF FAST %K

    # FUNCTION:

    # FPD:
    TS = is.timeSeries(close)
    fpd = emaTA(fpkTA(close, high, low, lag1), lag2)
    if (TS) colnames(fpd)<-"FPD"

    # Return Value:
    fpd
}


# ------------------------------------------------------------------------------


spdTA =
function(close, high, low, lag1 = 5, lag2 = 3, lag3 = 9)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Return Slow %D Indicator

    # Details:
    #   SLOW %D INDICATOR:
    #   EMA OF FAST %D

    # FUNCTION:

    # SPD:
    TS = is.timeSeries(close)
    spd = emaTA(fpdTA(close, high, low, lag1, lag2), lag3)
    if (TS) colnames(spd)<-"SPD"

    # Return Value:
    spd
}


# ------------------------------------------------------------------------------


apdTA =
function(close, high, low, lag1 = 5, lag2 = 3, lag3 = 9, lag4 = 9)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns Averaged %D Indicator

    # Details:
    #   AVERAGED %D INDICATOR: EMA OF SLOW %D

    # FUNCTION:

    # APD:
    TS = is.timeSeries(close)
    apd = emaTA(spdTA(close, high, low, lag1, lag2, lag3), lag4)
    if (TS) colnames(apd)<-"APD"

    # Return Value:
    apd
}


# ------------------------------------------------------------------------------


wprTA =
function(close, high, low, lag = 50)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns Williams %R Indicator

    # Details:
    #   Short Term: 5 to 49 days
    #   Intermediate Term: 50 to 100 day
    #   Long Term: More than 100 days

    # FUNCTION:

    # Check:
    TS = is.timeSeries(close)
    if (TS) {
        X = close
        close = as.vector(close)
        high = as.vector(high)
        low = as.vector(low)
    }

    # %R:
    minlag =
    function(x, lag)
    {
        xm = x
        for (i in 1:lag){
            x1 = c(x[1], x[1:(length(x)-1)])
            xm = pmin(xm, x1)
            x = x1
        }
        xm
    }
    maxlag =
    function(x, lag)
    {
        xm = x
        for (i in 1:lag){
            x1 = c(x[1], x[1:(length(x)-1)])
            xm = pmax(xm, x1)
            x = x1
        }
        xm
    }
    xmin = minlag(low, lag)
    xmax = maxlag(high, lag)
    wpr = (close - xmin) / (xmax -xmin)

    # to timeSeries:
    if (TS) {
        wpr = matrix(wpr)
        rownames(wpr) = rownames(series(X))
        colnames(wpr) = "WPR"
        series(X) = wpr
    } else {
        X = wpr
    }


    # Return Result:
    X
}


# ------------------------------------------------------------------------------


rsiTA =
function(close, lag = 14)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns Relative Strength Index Indicator

    # FUNCTION:

    # Check:
    TS = is.timeSeries(close)
    if (TS) {
        X = close
        close = as.vector(close)
    }

    # RSI:
    sumlag =
    function(x, lag){
        xs = x
        for (i in 1:lag){
            x1 = c(x[1],x[1:(length(x)-1)])
            xs = xs + x1
            x = x1
        }
        xs
    }
    close1 = c(close[1],close[1:(length(close)-1)])
    x = abs(close - close1)
    x[close<close1] = 0
    rsi = sumlag(x, lag)/sumlag (abs(close-close1), lag)
    rsi[1] = rsi[2]

    # to timeSeries:
    if (TS) {
        rsi = matrix(rsi)
        rownames(rsi) = rownames(series(X))
        colnames(rsi) = "RSI"
        series(X) = rsi
    } else {
        X = rsi
    }

    # Return Result:
    X
}


################################################################################
# MORE INDICATORS:


accelTA =
function(x, n = 12)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical Indicator Acceleration

    # FUNCTION:

    # Check:
    TS = is.timeSeries(x)
    if (TS) {
        stopifnot(isUnivariate(x))
        X = x
        x = as.vector(x)
    }

    # Indicator:
    accel = diff( x[(n+1):length(x)]-x[1:(length(x)-n)] )
    accel = c(rep(NA, n+1), accel)
    if (TS) {
        accel = matrix(accel)
        rownames(accel) = rownames(series(X))
        colnames(accel) = "ACCEL"
        series(X) = accel
    } else {
        X = accel
    }

    # Return Value:
    X
}


# ------------------------------------------------------------------------------


adiTA =
function(high, low, close, volume)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical indicator accumulation/distribution

    # FUNCTION:

    # Indicator:
    adi = (2 * close - high - low) / (high - low) * volume
    ADI = cumsum(as.vector(adi))

    # Time Series ?
    if (is.timeSeries(adi)) {
        ADI = matrix(ADI)
        colnames(ADI) = "ADI"
        rownames(ADI) = rownames(series(high))
        X = high
        series(X) = ADI
    } else {
        X = ADI
    }

    # Return Value:
    X
}


# ------------------------------------------------------------------------------


adoscillatorTA =
function(open, high, low, close)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical indicator Accumulation/Distribution Oscillator

    # FUNCTION:

    # Indicator:
    X = (high - open + close - low) / (high - low) * 50

    # Time Series ?
    if (is.timeSeries(open)) colnames(X)<-"ADO"

    # Return Value:
    X
}


# ------------------------------------------------------------------------------


bollingerTA =
function(x, lag = 20, n.sd = 2)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical Indicator Bollinger Band

    # FUNCTION:

    # Check:
    TS = is.timeSeries(x)
    if (TS) {
        X = x
        stopifnot(isUnivariate(x))
        x = as.vector(x)
    }

    # Indicator:
    n = lag
    mean = c(rep(NA, n-1), SMA(x = x, n = n))
    std = c(rep(NA, n-1), n.sd*sqrt(rollVar(x = x, n = n)))
    bollinger = cbind(UPPER = mean+std, BOLLINGER = x, LOWER = mean-std)

    if (TS) {
        rownames(bollinger) = rownames(series(X))
        series(X) = bollinger
    } else {
        rownames(bollinger) = as.character(1:length(x))
        X = bollinger
    }

    # Return Value:
    X
}


# ------------------------------------------------------------------------------


chaikinoTA =
function(high, low, close, volume, lag1 = 10, lag2 = 3)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical indicator Chaikin's oscillator

    # FUNCTION:

    # Indicator:
    adi = adiTA(high, low, close, volume)

    chaikino =
        emaTA(adi, lambda = 2/(lag1+1), startup = 0) -
        emaTA(adi, lambda = 2/(lag2+1), startup = 0 )

    if(is.timeSeries(chaikino)) colnames(chaikino)<-"CHAIKINO"

    # Return Value:
    chaikino
}


# ------------------------------------------------------------------------------


chaikinvTA =
function(high, low, lag1 = 10, lag2 = 10)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical indicator Chaikin's volatility

    # Arguments:
    #   lag1 - length of EMA range
    #   lag2 - length of change

    # FUNCTION:

    # Indicator:
    RT = emaTA(high-low, 2/(lag1+1), startup = 0)
    if (is.timeSeries(RT)) rt = as.vector(RT)
    chaikinv = (rt[-(1:lag2)]/rt[1:(length(rt)-lag2)]-1)*100
    chaikinv = c(rep(NA, lag2), chaikinv)

    if (is.timeSeries(RT)) {
        x = matrix(chaikinv)
        rownames(x) = rownames(series(RT))
        colnames(x) = "CHAIKINV"
    } else {
        x = chaikinv
    }

    # Return Value:
    x
}


# ------------------------------------------------------------------------------


garmanklassTA =
function(open, high, low, close)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the Garman-Klass estimate of volatility

    # FUNCTION:

    # Check:
    TS = is.timeSeries(open)
    if (TS) {
        x = open
        open = as.vector(open)
        high = as.vector(high)
        low = as.vector(low)
        close = as.vector(close)
    }

    # Indicator:
    prices = log(cbind(open, high, low, close))
    n = nrow(prices)
    alpha = 0.12
    f = 0.192
    u = high-open
    d = low-open
    cc = close - open
    oc = (prices[2:n, 1] - prices[1:(n - 1), 4])^2
    garmanklass = 0.511*(u-d)^2 - 0.019*(cc*(u+d) - 2*u*d) - 0.383*cc^2
    garmanklass = sqrt(((1 - alpha)*garmanklass[2:n])/(1-f) + (alpha*oc)/f)
    garmanklass = c(NA, garmanklass)

    if (TS) {
        garmanklass = matrix(garmanklass)
        colnames(garmanklass) = "GK"
        rownames(garmanklass) = rownames(series(x))
        series(x) = garmanklass
    } else {
        x = garmanklass
    }

    # Return Value:
    x
}


# ------------------------------------------------------------------------------


nviTA =
function(close, volume)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Reurns the technical indicator negative volume index

    # FUNCTION:

    # Check:
    TS = is.timeSeries(close)
    if (TS) {
        x = close
        close = as.vector(close)
        volume = as.vector(volume)
    }

    # Indicator:
    ind = rep(0, length(close)-1)
    ind[diff(volume) < 0] = 1
    ch = rocTA(close, lag = 1) / 100
    nvi = cumsum(ch * c(0, ind))

    # Time Series Output ?
    if (TS) {
        nvi = matrix(nvi)
        colnames(nvi) = "NVI"
        rownames(nvi) = rownames(series(x))
        series(x) = nvi
    } else {
        x = nvi
    }

    # Return Value:
    x
}



# ------------------------------------------------------------------------------


obvTA =
function(close, volume)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical indicator on balance volume

    # FUNCTION:

    # Check:
    TS = is.timeSeries(close)
    if (TS) {
        x = close
        close = as.vector(close)
        volume = as.vector(volume)
    }

    # Indicator:
    obv = cumsum(volume * c(0, sign(diff(close))))

    # Time Series Output ?
    if (TS) {
        obv = matrix(obv)
        colnames(obv) = "OBV"
        rownames(obv) = rownames(series(x))
        series(x) = obv
    } else {
        x = obv
    }

    # Return Value:
    x
}



# ------------------------------------------------------------------------------


pviTA =
function(close, volume)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical indicator positive volume

    # FUNCTION:

    # Check:
    TS = is.timeSeries(close)
    if (TS) {
        x = close
        close = as.vector(close)
        volume = as.vector(volume)
    }

    # Indicator:
    ind = rep(0, length(close)-1)
    ind[diff(volume) > 0] = 1
    ch = rocTA(close, lag = 1)/100
    pvi = cumsum(ch * c(0, ind))

    # Time Series Output ?
    if (TS) {
        pvi = matrix(pvi)
        colnames(pvi) = "PVI"
        rownames(pvi) = rownames(series(x))
        series(x) = pvi
    } else {
        x = pvi
    }

    # Return Value:
    x
}


# ------------------------------------------------------------------------------


pvtrendTA =
function(close, volume)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical indicator price-volume trend

    # FUNCTION:

    # Check:
    TS = is.timeSeries(close)
    if (TS) {
        x = close
        close = as.vector(close)
        volume = as.vector(volume)
    }

    # Indicator:
    m = length(close)
    ch = cumsum( volume * c(0, (close[2:m]/close[1:(m-1)]-1)*100))

    # Time Series Output ?
    if (TS) {
        ch = matrix(ch)
        colnames(ch) = "PVTREND"
        rownames(ch) = rownames(series(x))
        series(x) = ch
    } else {
        x = ch
    }

    # Return Value:
    x
}


# ------------------------------------------------------------------------------


williamsadTA =
function(high, low, close)
{   # A function written by Diethelm Wuertz

    # Description:
    #

    # FUNCTION:

    # Check:
    TS = is.timeSeries(high)
    if (TS) {
        x = high
        high = as.vector(high)
        low = as.vector(low)
        close = as.vector(close)
    }

    # Indicator:
    ind = c(0, sign(diff(close)))
    williamsad = vector("numeric", length(close))
    ind.pos = (ind == 1)
    ind.neg = (ind == -1)
    williamsad[ind.pos] = (close - low)[ind.pos]
    williamsad[ind.neg] =  - (high - close)[ind.neg]
    williamsad = cumsum(williamsad)
    names(williamsad) = as.character(1:length(x))

    # Time Series Output ?
    if (TS) {
        williamsad = matrix(williamsad)
        colnames(williamsad) = "WAD"
        rownames(williamsad) = rownames(series(x))
        series(x) = williamsad
    } else {
        x = williamsad
    }

    # Return Value:
    x
}



# ------------------------------------------------------------------------------


williamsrTA =
function(high, low, close, lag = 20)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Returns the technical indicator Willimas' accumulation/distribution

    # FUNCTION:

    # Check:
    TS = is.timeSeries(high)
    if (TS) {
        x = high
        high = as.vector(high)
        low = as.vector(low)
        close = as.vector(close)
    }

    # Indicator:
    n = lag
    hh = rollFun(high, n, trim = FALSE, FUN = max)
    ll = rollFun(low, n, trim = FALSE, FUN = min)
    williamsr = (hh-close)/(hh-ll)*100
    names(williamsr) = as.character(1:length(x))

    # Time Series Output ?
    if (TS) {
        williamsr = matrix(williamsr)
        colnames(williamsr) = "WR"
        rownames(williamsr) = rownames(series(x))
        series(x) = williamsr
    } else {
        x = williamsr
    }

    # Return Value:
    williamsr
}


################################################################################


SMA =
function(x, n = 5)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes a Simple Moving Average

    # Arguments:
    #   x - an univariate object of class "timeSeries" or a numeric vector

    # FUNCTION:

    # Rolling Mean:
    ans = rollFun(x = x, n = n, FUN = mean)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


EWMA =
function(x, lambda, startup = 0)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes an Exponentially Weighted Moving Average

    # Arguments:
    #   x - an univariate object of class "timeSeries" or a numeric vector

    # FUNCTION:

    # EWMA:
    ans = emaTA(x = x, lambda = lambda, startup = startup)

    # Return Value:
    ans
}


################################################################################


.dailyTA =
function(X, indicator = "ema", select = "Close", lag = 9)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute an indicator for technical analysis.

    # Arguments:
    #   X - a data.frame or timeSeries object with named
    #       columns, at least "Close", "Open", "High" and
    #       "Low". The order of the columns may be arbitrary.

    # FUNCTION:

    if (is.timeSeries(X)) {
        x = series(X)
    } else {
        stop("X must be a timeSeries object!")
    }

    if (indicator == "ema") {
        ans = emaTA(x = x[, select], lambda = lag)
    }

    if (indicator == "bias") {
        ans = biasTA(x = x[, select], lag = lag)
    }

    if (indicator == "medprice") {
        ans = medpriceTA(high = x[, "High"], low = x[, "Low"])
    }

    if (indicator == "typicalprice") {
        ans = typicalpriceTA(high = x[, "High"], low = x[, "Low"],
            close = x[, "Close"])
    }

    if (indicator == "wclose") {
        ans = wcloseTA(high = x[, "High"], low = x[, "Low"],
            close = x[, "Close"])
    }

    if (indicator == "roc") {
        ans = rocTA(x = x[, select], lag = lag)
    }

    if (indicator == "osc") {
        if (length(lag) < 2)
            stop("At least two lags must be specified!")
        ans = oscTA(x = x[, select], lag1 = lag[1], lag2 = lag[2])
    }

    if (indicator == "mom") {
        ans = momTA(x = x[, select], lag = lag)
    }

    if (indicator == "macd") {
        if (length(lag) < 2)
            stop("At least two lags must be specified!")
        ans = macdTA(x = x[, select], lag1 = lag[1], lag2 = lag[2])
    }

    if (indicator == "cds") {
        if (length(lag) < 3)
            stop("At least three lags must be specified!")
        ans = cdsTA(x = x[, select], lag1 = lag[1], lag2 = lag[2],
            lag3 = lag[3])
    }

    if (indicator == "cdo") {
        if (length(lag) < 3)
            stop("At least three lags must be specified!")
        ans = cdoTA(x = x[, select], lag1 = lag[1], lag2 = lag[2],
            lag3 = lag[3])
    }

    if (indicator == "vohl") {
        ans = vohlTA(high = x[, "High"], low = x[, "Low"])
    }

    if (indicator == "vor") {
        ans = vorTA(high = x[, "High"], low = x[, "Low"])
    }

    if (indicator == "fpk") {
        ans = fpkTA(close = x[, "Close"], high = x[, "High"],
            low = x[, "Low"], lag = lag)
    }

    if (indicator == "fpd") {
        if (length(lag) < 2)
            stop("At least two lags must be specified!")
        ans = fpdTA(close = x[, "Close"], high = x[, "High"],
            low = x[, "Low"], lag1 = lag[1], lag2 = lag[2])
    }

    if (indicator == "spd") {
        if (length(lag) < 3)
            stop("At least three lags must be specified!")
        ans = spdTA(close = x[, "Close"], high = x[, "High"],
            low = x[, "Low"], lag1 = lag[1], lag2 = lag[2],
            lag3 = lag[3])
    }

    if (indicator == "apd") {
        if (length(lag) < 4)
            stop("At least four lags must be specified!")
        ans = apdTA(close = x[, "Close"], high = x[, "High"],
            low = x[, "Low"], lag1 = lag[1], lag2 = lag[2],
            lag3 = lag[3], lag4 = lag[4])
    }

    if (indicator == "wpr") {
        ans = wprTA(close = x[, "Close"], high = x[, "High"],
            low = x[, "Low"], lag = lag)
    }

    if (indicator == "rsi") {
        ans = rsiTA(close = x[, "Close"], lag = lag)
    }

    # Return Value:
    timeSeries(data = matrix(ans, ncol = 1), charvec = X@positions,
        units = indicator, format = "ISO", zone = "GMT",
        FinCenter = "GMT")
}


################################################################################


.tradeSignals =
function(Positions)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes trade signals from trading positions

    # FUNCTION:

    # Get Signals from Positions:
    stopifnot(is.timeSeries(Positions))
    Signals = diff(Positions, pad = 0)/2
    Signals = Signals[abs(series(Signals)) == 1,]

    # Return Value:
    Signals
}


# ------------------------------------------------------------------------------


.tradeLengths =
function(tradeSignals)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes trade length from trading signals

    # FUNCTION:

    # Get Lengths from Signals:
    stopifnot(is.timeSeries(tradeSignals))
    data = diff(time(tradeSignals))
    charvec = tradeSignals@positions[-1]
    tradeLengths = timeSeries(data, charvec, units = "tradeLengths")

    # Return Value:
    tradeLengths
}


# ------------------------------------------------------------------------------


.hitRate =
function(Returns, Positions)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes hit rates from returns and positions

    # FUNCTION:

    # Compute hit rate:
    Indicator = (Positions * sign(Returns) + 1) / 2
    Rate = mean ( as.vector(Indicator), na.rm = TRUE )
    Rate = round(Rate, digits = 3)

    # Return Value:
    Rate
}


################################################################################
