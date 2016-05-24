
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
# FUNCTION:                   DESCRIPTION:
#  backtestPlot                Creates a summary of backtesting plots
#   backtestAssetsPlot         Plots assets used in a portfolio backtest   
#   backtestWeightsPlot        Plots recommended weights from a backtest
#   backtestRebalancePlot      Plots rebalanced weights of a backtest 
#   backtestPortfolioPlot      Plots benchmark and portfolio series
#   backtestDrawdownPlot       Plots the drawdown of the portfolio backtest
#   backtestReportPlot         Prints backtest report
################################################################################


backtestPlot <-
    function(object, which = "all", labels = TRUE, legend = TRUE,
    at = NULL, format = NULL, cex=0.6, font=1, family="mono")
{
    # A function implemented by Diethelm Wuertz and William Chen

    # Description:
    #   Creates a summary of backtesting plots
    
    # Arguments:
    #   object - a list as returned by the function portfolioSmoothing()
    #   which - which plots should be displayed
    #   labels - a logical flag, should automated labels added to the plot
    
    # FUNCTION:
    
    # Frame:
    if (any(which == "all"))
    par(mfrow = c(3, 2), mar = c(1.5, 4, 5, 2), oma = c(5, 1, 0, 1))
       
    # Plot:
    if(any(which == "1") || which == "all")
        backtestAssetsPlot (object, labels, legend, at, format)
    if(any(which == "2") || which == "all")
        backtestWeightsPlot (object, labels, legend, at, format)
    if(any(which == "3") || which == "all")
       backtestRebalancePlot (object, labels, legend, at, format)
    if(any(which == "4") || which == "all")
        backtestPortfolioPlot(object, labels, legend, at, format)
    if(any(which == "5") || which == "all")
        backtestDrawdownPlot(object, labels, legend, at, format)
    if(any(which == "6" )|| which == "all")
        backtestReportPlot(object, cex=cex, font=font, family=family)
        
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------
# Plot 1:

   
backtestAssetsPlot <-
    function(object, labels=TRUE, legend=TRUE, at=NULL, format=NULL)
{
    # A function implemented by Diethelm Wuertz and William Chen
    
    # Description:
    #   Plots assets used in a portfolio backtest
    
    # Arguments:
    #   object - a list as returned by the function portfolioSmoothing()
    #   labels - a logical flag, should automated labels added to the plot
    
    # FUNCTION:
    
    # Settings:
    data <- object$data
    benchmark <- object$benchmarkName
    assets <- object$assetsNames
    
    # Time Axis:
    if (is.null(at)) at <- paste(unique(atoms(time(data))[,1]), "12-31", sep="-")
    if (is.null(format)) Format <- "%b/%y" else Format <- format
    
    # Labels ?
    if (labels) {
        main <- "Index Series"
        xlab <- ""
        ylab <- "Cumulated log Returns"
    } else {
        main <- ""
        xlab <- ""
        ylab <- ""  
    }
    
    # Series:
    X <- data[, benchmark]
    
    # ylim - Plot Range:
    nAssets <- length(assets)
    MAX <- -1.0e99
    for (i in 1:nAssets) MAX = max(c(MAX, cumsum(data[, assets[i]])) )
    MAX <- max(MAX, cumsum(data[, benchmark]))
    MIN <- 1.0e99
    for (i in 1:nAssets) MIN <- min(MIN, cumsum(data[, assets[i]]))
    MIN <- min(MIN, cumsum(data[, benchmark]))
    rangeY <- c(MIN, MAX)
    
    # xlim - Plot Range:
    xlim <- range(time(colCumsums(data[, benchmark])))
    shift <- round(0.20 *as.integer(diff(xlim)), 0) * 24 * 60 * 60
    rangeX <- c(round(xlim[1]-shift), xlim[2])
    Days <- 1:as.integer(diff(xlim))
    Time <- as.character(xlim[1] + Days*24*60*60)
    range.tS <- timeSeries(data = matrix(rep(0, length(Time))), as.character(Time))
    
    # Limits:
    xlim <- rangeX
    ylim <- rangeY
    
    # Plot:
    plot(X, type = "n", xaxt = "n", at = at, format = Format,
        xlim = xlim, ylim = ylim, main = "", xlab = "", ylab = "") 
    grid(NA, ny = NULL)
    abline(v = as.POSIXct(at), lty = 3, col = "brown")
    
    # Add Lines:
    lines(colCumsums(data[, benchmark]), col = "black")
    lines(colCumsums(data[, benchmark]), col = "black")
    for (i in 1:nAssets) lines( colCumsums(data[, assets[i]]), col = i+1)
    
    # Asset Names:
    Benchmark <- abbreviate(benchmark, 4)
    Assets <- abbreviate(assets, 4)
    assetsList <- c(Benchmark, Assets)
    assetsTitle <- paste(Benchmark, " ~ ", 
        paste( Assets, collapse = " - ", sep = ""), sep="")
    
    # Add Title:
    if (labels) {
        title(main = main, xlab = xlab, ylab = ylab)
    }
    
    # Add Legend and Subtitle:
    if (legend) {
        mtext(assetsTitle, line = 0.5, cex = 0.7)
        legend("topleft",
            legend = assetsList,
            bty = "n",
            text.col = 1:(nAssets+1),
            cex = 0.8)
    }
   
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------
# Plot 2:


backtestWeightsPlot <-
    function(object, labels=TRUE, legend=TRUE, at=NULL, format=NULL)
{
    # A function implemented by Diethelm Wuertz and William Chen
    
    # Description:
    #   Plots recommended weights from a portfolio backtest
   
    # Arguments:
    #   object - a list as returned by the function portfolioSmoothing()
    #   labels - a logical flag, should automated labels added to the plot

    # FUNCTION:

    # Settings:
    data <- object$data
    weights <- object$smoothWeights
    assets <- object$assetsNames
    benchmark <- object$benchmarkName
    horizon <- getWindowsHorizon(object$backtest)
    smoothing <- getSmootherLambda(object$backtest)
    startup <- "1m"
    horizonLength = as.numeric(substr(horizon, 1, nchar(horizon)-1))
    horizonUnit = substr(horizon, nchar(horizon), nchar(horizon))
    stopifnot(horizonUnit == "m")
    
    # Time Axis:
    if (is.null(at)) at <- paste(unique(atoms(time(data))[,1]), "12-31", sep="-")
    if (is.null(format)) Format <- "%b/%y" else Format <- format
    
    # Labels ?
    if (labels) {
        main <- "Weights Recommendation"
        xlab <- ""
        ylab <- "Asset Weights %"
    } else {
        main <- ""
        xlab <- ""
        ylab <- ""
    }

    # Series:
    X <- data[, benchmark]
    nAssets <- length(assets)
    naWeights <- matrix(rep(NA, times=horizonLength*nAssets), ncol=nAssets)

    # Lmits:
    xlim <- range(data)
    ylim <- c(0, 100)
    
    # Plot:
    plot(X, type = "n", xaxt = "n", las = 2, at = at, format = Format,
        xlim = xlim, ylim = ylim, main = "", xlab = "", ylab = "") 
    grid(NA, ny = NULL)
    abline(v = as.POSIXct(at), lty = 3, col = "brown")
    
    # Add Lines:
    lines(X, col = "black")
    tS <- 100 * timeSeries(weights)
    for (i in 1:nAssets) lines(tS[, i], col = i+1)

    # Asset Names:
    Benchmark <- abbreviate(benchmark, 4)
    Assets <- abbreviate(assets, 4)
    assetsList <- c(Benchmark, Assets)
    assetsTitle <- paste(Benchmark, " ~ ", 
        paste(Assets, collapse = " - ", sep = ""), sep="")
    
    # Add Title:
    if (labels) {
        title(main = main, xlab = xlab, ylab = ylab)
        text <- paste(
            "Horizon = ", horizon,
            "| Smoothing:", smoothing,
            "| Startup:", startup,
            "| Shift 1m")
        mtext(text, line = 0.5, cex = 0.7)
    }
    
    # Add Legend:
    if (legend) {
        legend("topleft",
            legend = assetsList,
            bty = "n",
            text.col = 1:(nAssets+1),
            cex = 0.8)
    }

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------
# Plot 3:


backtestRebalancePlot <-
    function(object, labels=TRUE, legend=TRUE, at=NULL, format=NULL)
{
    # A function implemented by Diethelm Wuertz and William Chen
    
    # Description:
    #   Plots rebalanced weights of a backtest

    # Arguments:
    #   object - a list as returned by the function portfolioSmoothing()
    #   labels - a logical flag, should automated labels added to the plot
   
    # FUNCTION:

    # Settings:
    data <- object$data
    weights <- object$smoothWeights
    assets <- object$assetsNames
    benchmark <- object$benchmarkName
    horizon <- getWindowsHorizon(object$backtest)
    smoothing <- getSmootherLambda(object$backtest)
    startup <- "1m"
    horizonLength <- as.numeric(substr(horizon, 1, nchar(horizon)-1))
    horizonUnit <- substr(horizon, nchar(horizon), nchar(horizon))
    stopifnot(horizonUnit == "m")
    horizon <- horizonLength
    
    # Time Axis:
    if (is.null(at)) at <- paste(unique(atoms(time(data))[,1]), "12-31", sep="-")
    if (is.null(format)) Format <- "%b/%y" else Format <- format
    
    # Labels ?
    if (labels) {
        main <- "Weights Rebalance"
        xlab <- ""
        ylab <- "Weights Changes %"
    } else {
        main <- ""
        xlab <- ""
        ylab <- ""
    }
   
    # Series:
    X <- data[, benchmark]
    nAssets <- length(assets)
    naWeights <- matrix(rep(NA, times = horizon * nAssets), ncol = nAssets)
    naWeights <- rbind(naWeights, rep(NA, times = nAssets))
    diffWeights <- rbind(naWeights, diff(weights))
    absSum <- function(x) { sum(abs(x)) }
    diffWeights <- apply(diffWeights, 1, FUN = absSum)
    diffWeights <- cbind(diffWeights, rbind(naWeights, diff(weights)))
    tS <- 100 * timeSeries(diffWeights[-seq(horizon + 1),],
        charvec = rownames(diffWeights)[-seq(horizon + 1)])   
        
    # Limits:
    xlim <- range(time(data))
    ylim <- range(tS)
    
    # Plot:
    plot(X, type = "n", xaxt = "n", las = 2, at = at, format = Format,
        xlim = xlim, ylim = ylim, main = "", xlab = "", ylab = "") 
    grid(NA, ny = NULL)
    abline(v = as.POSIXct(at), lty = 3, col = "brown")
    abline(h=0, col="darkgrey")
    
    # Add Lines:
    # lines(X)
    lines(tS[, 1], type = "h", lwd = 1, col = "darkgrey")
    for (i in 2:NCOL(tS)) lines(tS[, i], col = i)
   
    # Asset Names:
    Benchmark <- abbreviate(benchmark, 4)
    Assets <- abbreviate(assets, 4)
    assetsList <- c(Benchmark, Assets)
    assetsTitle <- paste(Benchmark, " ~ ", 
        paste(Assets, collapse = " - ", sep = ""), sep="")
    
    # Add Title:
    if (labels) {
        title(main = main, xlab = xlab, ylab = ylab)
        text <- paste(
            "Horizon = ", horizon,
            "| Smoothing:", smoothing,
            "| Startup:", startup,
            "| Shift 1m")
        mtext(text, line = 0.5, cex = 0.7)
        # mText = paste("Start:", rownames(object$smoothWeights)[1])
        # mtext(mText, side = 4, line = 0, adj = 0, col = "darkgrey", cex = 0.65)
    }
    
    # Add Legend:
    if (legend) {
        legend("topleft",
            legend = assetsList,
            bty = "n",
            text.col = 1:(nAssets+1),
            cex = 0.8)
    }

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------
# Plot 4:


backtestPortfolioPlot <-
    function(object, labels=TRUE, legend=TRUE, at=NULL, format=NULL)
{
    # A function implemented by Diethelm Wuertz and William Chen
    
    # Description:
    #   Plots daily, benchmark and portfolio series of a portfolio backtest

    # Arguments:
    #   object - a list as returned by the function portfolioSmoothing()
    #   labels - a logical flag, should automated labels added to the plot

    # FUNCTION:

    # Settings:
    data <- object$data
    portfolioReturns <- object$portfolioReturns
    benchmarkReturns <- object$benchmarkReturns
    assets <- object$assetsNames
    benchmark <- object$benchmarkName
    horizon <- getWindowsHorizon(object$backtest)
    smoothing <- getSmootherLambda(object$backtest)
    startup <- "1m"
    offsetReturn <- object$offsetReturn
    
    # Time Axis:
    if (is.null(at)) at <- paste(unique(atoms(time(data))[,1]), "12-31", sep="-")
    if (is.null(format)) Format <- "%b/%y" else Format <- format
    
    # Labels ?
    if (labels) {
        main <- "Portfolio vs Benchmark"
        xlab <- ""
        ylab <- "Cumulated log Returns"   
    } else {
        main <- ""
        xlab <- ""
        ylab <- ""  
    }
   
    # Series:
    X <- data[, benchmark]

    # Cumulated Return Series:
    cumX <- colCumsums(X)
    cumP <- portfolioReturns + offsetReturn
    cumB <- benchmarkReturns + offsetReturn
    offsetTS <- timeSeries(offsetReturn, charvec = names(offsetReturn),
        units = "offsetReturn")
    cumP <- rbind(offsetTS, cumP)
    cumB <- rbind(offsetTS, cumB)
    MAX <- max(as.vector(series(cumP)), as.vector(series(cumB)),
        as.vector(series(cumX)))
    MIN <- min(as.vector(series(cumP)), as.vector(series(cumB)),
        as.vector(series(cumX)))
       
    # Limits:
    xlim <- c(as.POSIXct(start(X)), as.POSIXct(end(X)))
    ylim <- c(MIN, MAX)
     
    # Plot:
    plot(X, type = "n", xaxt = "n", at = at, format = Format,
        xlim = xlim, ylim = ylim, main = "", xlab = "", ylab = "") 
    grid(NA, ny = NULL)
    abline(v = as.POSIXct(at), lty = 3, col = "brown")
    abline(h=0, col="darkgrey")
    
    # Add Lines:
    lines(cumX, col = "black")
    lines(cumP-cumB, type = "h", col = "grey")
    lines(cumP, col = "red", lwd = 2)
    lines(cumB, col = "blue", lwd = 2)

    # Asset Names:
    Benchmark <- abbreviate(benchmark, 4)
    Assets <- abbreviate(assets, 4)
    assetsList <- c(Benchmark, Assets)
    assetsTitle <- paste(Benchmark, " ~ ", 
        paste(Assets, collapse = " - ", sep = ""), sep="")
    nAssets <- length(assetsList)
    
    # Add Title:
    if (labels) {
        title(main = main, xlab = xlab, ylab = ylab)
        text <- paste(
            "Horizon = ", horizon,
            "| Smoothing:", smoothing,
            "| Startup:", startup,
            "| Shift 1m")
        mtext(text, line = 0.5, cex = 0.7)
        # mText <- Type = getType(object$spec)
        # Estimator <- getEstimator(object$spec)
        # if (Type == "MV") mText = paste(mText, "|", Estimator)
        # mtext(mText, side = 4, line = 0, adj = 0, col = "darkgrey", cex = 0.7)
    }
    
    # Add Legend:
    if (legend) {
        legend("topleft",
            legend = assetsList,
            bty = "n",
            text.col = 1:(nAssets+1),
            cex = 0.8)
    }
    
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------
# Plot 5:


backtestDrawdownPlot <- 
    function(object, labels=TRUE, legend=TRUE, at=NULL, format=NULL)
{
    # A function implemented by Diethelm Wuertz and William Chen

    # Description:
    #   Creates Backtest Portfolio Plot

    # Arguments:
    #   object - a list as returned by the function portfolioSmoothing()
    #   labels - a logical flag, should automated labels added to the plot

    # FUNCTION:
    
    # Settings:
    data <- object$data
    assets <- object$assetsNames
    benchmark <- object$benchmarkName
    horizon <- getWindowsHorizon(object$backtest)
    smoothing <- getSmootherLambda(object$backtest)
    startup <- getSmootherInitialWeights(object$backtest)
    weights <- as.timeSeries(object$smoothWeights)
    
    # Align Data:
    Data <- .align.timeSeries(data)/100
    
    # Time Axis:
    if (is.null(at)) at <- paste(unique(atoms(time(data))[,1]), "12-31", sep="-")
    if (is.null(format)) Format <- "%b/%y" else Format <- format
    
    # Labels ?
    if (labels) {
        main <- "Portfolio vs Benchmark"
        xlab <- ""
        ylab <- "Drawdown"
    } else {
        main <- ""
        xlab <- ""
        ylab <- ""
    }
    
    # Extract the Time Stamps:
    tS <- time(Data)
    tW <- time(weights)
        
    # Problem when rebalance day lands on a Weekend - 
    #   need to change the date to the nearest Monday
    if (any(isWeekend(tW))){
        weekend.tW <- tW[isWeekend(tW)]       
        # WC: check timeNdayOnOrAfter function, the nday = 2 is a Monday!???
        tW <- sort(c(tW[!isWeekend(tW)], timeNdayOnOrAfter(weekend.tW, 2)))
        # replace old times with new times
        time(weights) <- tW
    }
            
    # Extract the Updated Revalance Dates:  
    Dates <- time(weights)
    
    # Subsetting the Data:
    data <- window(Data, start(weights), end(weights))
    
    # Check whether we have data past the last balance date
    # i.e. last balance date won't take place if we don't have the return series
    if (end(data) < end(weights)){ 
        n <- length(Dates)-1 
    } else {n = length(Dates)
        Dates <- c(Dates, end(data))
    }
    
    # Calculate the portfolio returns for the given weights:
    # assume we start investing the new weights on the rebalance date
    pf <- NULL
    a <- NULL
    for (i in 1:n){
        temp <- window(data, Dates[i], Dates[i+1])[,assets]
        nr <- nrow(temp)
        if (i != n) temp = temp[-nr,]
        a <- c(a, nrow(temp))
        pf <- c(pf, pfolioReturn(temp, as.numeric(weights[i,])))
    }
    
    # Drawdown Plot Settings:
    stopifnot(length(pf) == length(rownames(data)))
    pf <- timeSeries(pf, charvec = rownames(data))
    pf.DD <- drawdowns(pf)
    benchmark.DD <- drawdowns(data[,benchmark]) 
    
    # Series:
    X <- Data[, benchmark]
    
    # Limits:
    xlim <- c(as.POSIXct(start(X)), as.POSIXct(end(X)))
    ylim <- range(c(pf.DD, benchmark.DD))
    
    # Plot:
    plot(X, type = "n", xaxt = "n", at = at, format = Format,
        xlim = xlim, ylim = ylim, main = "", xlab = "", ylab = "") 
    grid(NA, ny = NULL)
    abline(v = as.POSIXct(at), lty = 3, col = "brown")
    
    # Add Lines:
    lines(benchmark.DD,  col = "blue", lwd = 2)
    lines(pf.DD, col = "red", lwd = 2)
    
    # Asset Names:
    Benchmark <- abbreviate(benchmark, 4)
    Assets <- abbreviate(assets, 4)
    assetsList <- c(Benchmark, Assets)
    assetsTitle <- paste(Benchmark, " ~ ", 
        paste(Assets, collapse = " - ", sep = ""), sep="")
    
    # Add Title:
    if (labels) {
        title(main = main, xlab = xlab, ylab = ylab)
        text <- paste("(Max)", "Portfolio DD =", round(min(pf.DD),2),
            "|", "Benchmark DD =", round(min(benchmark.DD),2))
        mtext(text, line = 0.5, cex = 0.7)
    }
    
    # Add Legend:
    if (legend) {
        legend("bottomleft",
            legend = c("Benchmark", "Portfolio"),
            bty = "n",
            text.col = c("blue", "red"),
            cex = 0.8)
    }
       
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------
# Plot 6 - Report:


backtestReportPlot <-
    function(object, cex=0.6, font=1, family="mono")
{    
    # A function implemented by Diethelm Wuertz and William Chen
    
    # Description:
    #   Prints backtest report as graphical plot
    
    # Arguments:
    #   object - a list as returned by the function portfolioSmoothing()
    
    # FUNCTION:
    
    # Settings:
    CEX <- cex
    
    
    # Start Plot:
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
      
    # Vertical Adjustment:
    z <- -2
        
    TEXT <- paste("Strategy:", getStrategyFun(object$backtest))
    mtext(TEXT, side = 3, line =  z + 3, adj = 0, font=2, family="sans", 
          cex=CEX)
        
    TEXT <-  capture.output(round(object$stats, 2))
    mtext(TEXT[1], side = 3, line = z + +2, adj = 0, font=font, family="mono", cex=CEX)
    mtext(TEXT[2], side = 3, line = z + +1, adj = 0, font=font, family="mono", cex=CEX)
    mtext(TEXT[3], side = 3, line = z + +0, adj = 0, font=font, family="mono", cex=CEX) 
    mtext(TEXT[4], side = 3, line = z + -1, adj = 0, font=font, family="mono", cex=CEX)
    mtext(TEXT[5], side = 3, line = z + -2, adj = 0, font=font, family="mono", cex=CEX) 
        
    TEXT <- capture.output(object$spec)[c(2,3,4,5,8)]
    mtext("Portfolio Specification:", side = 3, line = z + -4, adj = 0, font=2, family="sans", cex=CEX)
        
    if (length(grep("CVaR",TEXT[2]))!=0) TEXT[2] = 
        gsub("CVaR", paste("CVaR |", getAlpha(object$spec)), TEXT[2])
       
    mtext(TEXT[2], side = 3, line = z + -5, adj = 0, font=font, family=family, cex=CEX)
    mtext(TEXT[3], side = 3, line = z + -6, adj = 0, font=font, family=family, cex=CEX)
    mtext(TEXT[4], side = 3, line = z + -7, adj = 0, font=font, family=family, cex=CEX)
    #text(TEXT[5], side = 3, line = z + -8, adj = 0, font=font, family=family, cex=CEX)
    
    TEXT <- capture.output(object$constraints)[1]      
    mtext("Constraints:", side = 3, line = z + -9, adj = 0, font=2, family="sans", cex=CEX)
    TEXT <- substr(TEXT[1], 4, 99)   
    mtext(TEXT, side = 3, line = z + -10, adj = 0, font=font, family=family, cex=CEX)
                        
    # Return Value:
    invisible()
 }


################################################################################

