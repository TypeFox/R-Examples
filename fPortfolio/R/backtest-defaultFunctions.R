
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
# FUNCTION:                DESCRIPTION:
#  equidistWindows          Defines default equal distant rolling windows
#  tangencyStrategy         Defines default tangency strategy portfolio
#  emaSmoother              Defines default EMA weights smoother
################################################################################


equidistWindows <-
    function(data, backtest = portfolioBacktest())
{
    # A function implemented by Diethelm Wuertz and William Chen
    
    # Description:
    #   Defines default equidistant rolling windows
   
    # Arguments:
    #   data - portfolio assets set, an object of class 'timeSeries'
    #   backtest - an object of class 'fPFOLIOBACKTEST'
   
    # Note:
    #   This is an example for a user defined windows function ...
   
    # Example:
    #   equidistWindows(as.timeSeries(data(LPP2005REC)))
   
    # FUNCTION:
   
    # Settings:
    horizon = getWindowsHorizon(backtest)
   
    # Rolling Windows:
    ans = rollingWindows(x = data, period = horizon, by = "1m")
   
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


tangencyStrategy <-
    function(data, spec = portfolioSpec(), constraints = "LongOnly", 
    backtest = portfolioBacktest())
{ 
    # A function implemented by Diethelm Wuertz and William Chen
    
    # FUNCTION:
    
    # Strategy Portfolio:
    strategyPortfolio <- try(tangencyPortfolio(data, spec, constraints))
    
    # If tangency portfolio doesn't exist take the minimum variance portfolio:
    if (class(strategyPortfolio) == "try-error") {
        strategyPortfolio <- minvariancePortfolio(data, spec, constraints)
    }
    
    # Return Value:
    strategyPortfolio
}


# ------------------------------------------------------------------------------


emaSmoother <-
    function(weights, spec, backtest)
{
    # A function implemented by Diethelm Wuertz and William Chen
    
    # Description:
    #   A user defined weights smoother for portfolio backtesting
   
    # Arguments:
    #   weights - a numeric matrix of weights
    #   spec - portfolio spec, an object of class fPFLOLIOSPEC
    #   backtest - portfolio backtest, an object of class fPFLOLIOBACKTEST
   
    # Example:
    #   ans = portfolioBacktesting( ... )
    #   emaSmoother(ans$weights, spec, backtest)
   
    # FUNCTION:
   
    # EMA Function:
    ema <- function (x, lambda) {
        x = as.vector(x)
        lambda = 2/(lambda + 1)
        xlam = x * lambda
        xlam[1] = x[1]
        ema = filter(xlam, filter = (1 - lambda), method = "rec")
        ema[is.na(ema)] <- 0
        as.numeric(ema) }
       
    # Lambda:
    lambda <- getSmootherLambda(backtest)
    lambdaLength <- as.numeric(substr(lambda, 1, nchar(lambda) - 1))
    lambdaUnit <- substr(lambda, nchar(lambda), nchar(lambda))
    stopifnot(lambdaUnit == "m")
    lambda <- lambdaLength
   
    # Initial Weights:
    nAssets <- ncol(weights)
    initialWeights = getSmootherInitialWeights(backtest)
    if (!is.null(initialWeights)) weights[1, ] = initialWeights

    # Compute Exponentially Smoothed Weights:
    smoothWeights1 = NULL
    for (i in 1:nAssets) {
        # print("first smooth")
        EMA <- ema(weights[, i], lambda = lambda)
        smoothWeights1 <- cbind(smoothWeights1, EMA)
    }
   
    # Double Smoothing ?
    doubleSmooth <- getSmootherDoubleSmoothing(backtest)
    if (doubleSmooth) {
        # print("second smooth")
        smoothWeights = NULL
        for (i in 1:nAssets) {
            EMA <- ema(smoothWeights1[, i], lambda = lambda)
            smoothWeights = cbind(smoothWeights, EMA)
        }
    } else {
        smoothWeights <- smoothWeights1
    }
   
    # Add Names:
    rownames(smoothWeights) <- rownames(weights)
    colnames(smoothWeights) <- colnames(weights)
   
    # Return Value:
    smoothWeights
}


################################################################################

