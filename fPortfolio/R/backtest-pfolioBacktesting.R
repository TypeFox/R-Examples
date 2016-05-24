
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


################################################################################
# FUNCTION:                    DESCRIPTION:
#  portfolioBacktesting         Performs a portfolio backtesting
#  portfolioSmoothing           Smoothes the weights of a portfolio backtesting
################################################################################


portfolioBacktesting <-
  function(formula, data, spec = portfolioSpec(), constraints = "LongOnly", 
           backtest = portfolioBacktest(), trace = TRUE)
  {
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Backtests a portfolio on rolling windows
    
    # Arguments:
    #   formula - a formula expression to select benchmark and assets
    #       from the data set
    #   data - data set of assets returns, an object of class fPFLOLIODATA 
    #       or timeSeries
    #   spec - portfolio specification, an object of class fPFLOLIOSPEC,
    #       by default as returned by the function portfolioSpec()
    #   constraints - portfolio constraints, a vector of character strings
    #   backtest - portfolio backtest specification, an object of 
    #       class fPFLOLIOBACKTEST, by default as returned by the function
    #       portfolioBacktest
    #   trace - a logical, should the backtesting be traced ?
    
    # Value:
    #   A list with the following elements 
    #   formula - the input formula 
    #   data - the input data set
    #   spec - the input portfolio specification 
    #   constraints - the input constraints 
    #   backtest - the input backtest specification  
    #   benchmarkName - the name of the benchmark returns 
    #   assetsNames - the names of the assets returns
    #   weights - the rolling weights matrix 
    #   strategyList - the rolling list of optimized portfolios 
    #   Sigma - ...
    
    # Details:
    #   Allows for user specified rolling Windows
    #   Smoothing is separated and can be user specified
    
    # Example:
    #   portfolioBacktesting(formula, data, spec, constraints, backtest)
    
    # FUNCTION:
    
    # Data:
    if (class(data) == "fPFOLIODATA") {
      Data <- data
      data <- getSeries(data)
    } else if (class(data) == "timeSeries") {
      Data <- portfolioData(data, spec)
    }
    
    # Constraints:
    if (class(constraints) == "fPFOLIOSPEC") {
      Constraints <- constraints
      constraints <- Constraints@stringConstraints
    } else if (class(constraints) == "character") {
      Constraints <- portfolioConstraints(data, spec, constraints)
    }
    
    # Formula, Benchmark and Asset Labels:
    benchmarkName = as.character(formula)[2]
    assetsNames <- strsplit(gsub(" ", "", as.character(formula)[3]), "\\+")[[1]]
    nAssets <- length(assetsNames)
    
    # Trace the Specifications and Data Info:
    if(trace) {
      cat("\nPortfolio Backtesting:\n")
      cat("\nAssets:             ", assetsNames)
      cat("\nBenchmark:          ", benchmarkName)
      cat("\nStart Series:       ", as.character(start(data)))
      cat("\nEnd Series:         ", as.character(end(data)))
      cat("\n  Type:             ", getType(spec))
      cat("\n  Cov Estimator:    ", getEstimator(spec))
      cat("\n  Solver:           ", getSolver(spec))
      cat("\nPortfolio Windows:  ", getWindowsFun(backtest))
      cat("\n  Horizon:          ", getWindowsHorizon(backtest))
      cat("\nPortfolio Strategy: ", getStrategyFun(backtest))
      cat("\nPortfolio Smoother: ", getSmootherFun(backtest))
      cat("\n  doubleSmoothing:  ", getSmootherDoubleSmoothing(backtest))
      cat("\n  Lambda:           ", getSmootherLambda(backtest))
    }
    
    # We invest in the "Strategy" or (return) efficient Portfolio:
    if(trace) {
      cat("\n\nPortfolio Optimization:")
      cat("\nOptimization Period\tTarget\tBenchmark\t Weights\n")
    }
    
    # Create Rolling Windows:
    windowsFun <- match.fun(getWindowsFun(backtest))
    rollingWindows <- windowsFun(data, backtest)
    from <- rollingWindows$from
    to <- rollingWindows$to
    
    # Roll the Portfolio:
    strategyFun <- match.fun(getStrategyFun(backtest))
    strategyList <- list()
    
    # WC: track the sigma over time:
    Sigma <- NULL
    
    for (i in 1:length(from)) 
    {
      # Optimize the Portfolio:
      pfSeries <- window(data[, assetsNames], start = from[i], end = to[i])
      bmSeries <- window(data[, benchmarkName], start = from[i], end = to[i])
      pfSeries <- portfolioData(pfSeries, spec)
      Sigma <- c(Sigma, mean(diag(getSigma(pfSeries))))
      strategy <- strategyFun(
        data = getSeries(pfSeries), 
        spec = spec, 
        constraints = constraints, 
        backtest = backtest)
      strategyList[[i]] <- strategy
      
      # Trace Optionally the Results:
      if (trace) {
        cat(as.character(from[i]), as.character(to[i]))
        spReturn <- getTargetReturn(strategy@portfolio)[[2]]
        cat("\t", round(spReturn[1], digits = 3))
        bmReturn <- mean(series(bmSeries))
        cat("\t", round(bmReturn, digits = 3))
        nAssets <- length(assetsNames)
        weights <- round(getWeights(strategy), digits = 3)
        cat("\t")
        for (i in 1:length(assetsNames)) cat("\t", weights[i])         
        cat("\t  * ", round(sum(weights), 2))
        cat("\n")
      }  
    }
    
    # Extract Portfolio Investment Weights for the current period:
    weights <- NULL
    for (i in 1:length(strategyList)) 
      weights <- rbind(weights, getWeights(strategyList[[i]]))
    rownames(weights) <- as.character(to)
    colnames(weights) <- assetsNames
    
    # Compose Result:
    ans <- list(
      formula = formula,
      data = data,
      spec = spec,
      constraints = constraints,
      backtest = backtest,
      benchmarkName = benchmarkName,
      assetsNames = assetsNames,
      weights = weights,
      strategyList = strategyList, 
      Sigma = Sigma)
    
    # Return Value:
    class(ans) <- c("portfolioBacktesting", "list")
    invisible(ans)
  }


# ------------------------------------------------------------------------------


portfolioSmoothing <-
  function(object, backtest=NULL, trace=TRUE)
  {
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Flexible Weights Smoother Function
    
    # Arguments:
    #   object - an object as returned by the function portfolioBacktesting()
    #   backtest - an S4 class object of 'FPFOLIOBACKTEST', the same as
    #       used in the function portfolioBacktesting() or a user modified
    #       version (obsolete)
    #   trace - a logical, should the computation be traced ?
    
    # Value:
    #   a list with the following entries
    
    # Example:
    #   data <- 100*SWX.RET; object=portfolioBacktesting(LP40~SBI+SPI+SII, data)
    #   portfolioSmoothing(object)
    
    # FUNCTION:
    
    # Obsolete Argument
    if (!is.null(backtest)) {
      warning("The backtest argument is obsolete and will be
            removed for the next release.")
    }
    
    # Backtest Settings:
    formula <- object$formula
    data <- object$data
    spec <- object$spec
    constraints <- object$constraints
    backtest <- object$backtest 
    benchmarkName <- object$benchmarkName
    assetsNames <- object$assetsNames
    weights <- object$weights
    skip <- getSmootherSkip(backtest)
    if (skip > 0) weights <- weights[-(1:skip), ]
    nAssets <- ncol(weights)
    
    # Add Smooth Weights to Backtest object:
    if (trace) print("smooth ...")
    smoother <- match.fun(getSmootherFun(backtest))
    smoothWeights <- object$smoothWeights <- smoother(weights, spec, backtest)
    
    # Compute Monthly Assets and Benchmark Returns:
    if (trace) print("aggregate ...")
    ow <- options("warn")
    options(warn = -1)
    monthlyAssets <- object$monthlyAssets <-
      applySeries(data[, assetsNames], by = "monthly", FUN = colSums)
    monthlyBenchmark <- object$monthlyBenchmark <-
      applySeries(data[, benchmarkName], by = "monthly", FUN = colSums)
    options(ow)   
    
    # Compute Offset Return of Rolling Portfolio compared to Benchmark:
    if (trace) print("offset ...")
    cumX <- colCumsums(data[, benchmarkName])
    lastX <- window(cumX, start = start(cumX), end = rownames(weights)[1] )
    offsetReturn <- as.vector(lastX[end(lastX), ])
    names(offsetReturn) <- as.character(end(lastX))
    object$offsetReturn <- offsetReturn
    
    # Backtest Return Series:
    Datum <- as.vector(rownames(smoothWeights))
    nDatum <- length(Datum)
    Portfolio = Benchmark = NULL
    for (i in 1:(nDatum-1)) {
      Portfolio <- rbind(Portfolio, as.vector((
        as.matrix(monthlyAssets[Datum[i+1], ]) %*% smoothWeights[Datum[i], ])))
      Benchmark <- rbind(Benchmark, as.vector(monthlyBenchmark[Datum[i+1], ]))
    }
    
    # Portfolio:
    P <- timeSeries(data = Portfolio, charvec = Datum[-1], units = "Portfolio")
    object$portfolioReturns <- portfolio <- colCumsums(P)
    object$P <- P
    
    # Benchmark:
    B <- timeSeries(data = Benchmark, charvec = Datum[-1], units = "Benchmark")
    object$benchmarkReturns <- benchmark <- colCumsums(B)
    object$B <- B
    
    daily <- colCumsums(data[, benchmarkName])
    Daily <- window(daily, start = start(portfolio), end = end(portfolio))
    
    portfolio <- portfolio - portfolio[1] + Daily[1]
    benchmark <- benchmark - benchmark[1] + Daily[1]
    
    # Add to backtest:
    object$portfolio <- portfolio
    object$benchmark <- benchmark
    
    # Backtest Statistics:
    P <- as.vector(P)
    B  <- as.vector(B)
    Stats <- c(sum(P, na.rm = TRUE), sum(B))
    Stats <- rbind(Stats, c(mean(P, na.rm = TRUE), mean(B)))
    Stats <- rbind(Stats, c(sd(P, na.rm = TRUE), sd(B)))
    Stats <- rbind(Stats, c(min(P, na.rm = TRUE), min(B)))
    colnames(Stats) <- c(
      "Portfolio",
      "Benchmark")
    rownames(Stats) <- c(
      "Total Return",
      "Mean Return",
      "StandardDev Return",
      "Maximum Loss")
    object$stats <- Stats
    
    # Return Value:
    class(object) <- c("portfolioSmoothing", "list")
    object
  } 


################################################################################

