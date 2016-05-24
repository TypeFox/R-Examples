
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
# FUNCTION:                    DESCRIPTION:
#  frontierPlot                 Plots efficient frontier
#   minvariancePoints             Adds minimum variance point
#   cmlPoints                     Adds market portfolio
#   cmlLines                      Adds capital market Line
#   tangencyPoints                Adds tangency portfolio point
#   tangencyLines                 Adds tangency line
#   equalWeightsPoints            Adds point of equal weights portfolio
#   singleAssetPoints             Adds points of single asset portfolios
#   twoAssetsLines                Adds EF for all combinations of two assets
#   sharpeRatioLines              Adds Sharpe ratio line
#   monteCarloPoints              Adds randomly produced feasible portfolios
# FUNCTION:                    DESCRIPTION:
#  frontierPlotControl          Sets frontier plot control parameters
# FUNCTION:                    DESCRIPTION:
#  tailoredFrontierPlot         Tailored frontier plot with addons
################################################################################


frontierPlot <-
function(object, frontier = c("both", "lower", "upper"),
    col = c("black", "grey"), add = FALSE, labels = TRUE,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, title = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots the efficient frontier
    
    # Arguments:
    
    # FUNCTION:

    # Check Settings:
    stopifnot(length(col) == 2)

    # Settings:
    frontier <- match.arg(frontier)
    fullFrontier = frontierPoints(object, frontier = "both",
        return = return, risk = risk, auto = auto)
    upperFrontier <- frontierPoints(object, frontier = "upper",
        return = return, risk = risk, auto = auto)
    lowerFrontier <- frontierPoints(object, frontier = "lower",
        return = return, risk = risk, auto = auto)

    # Check for 'xlim' Argument:
    Arg <- match.call(expand.dots = TRUE)
    m <- match(c("xlim", "ylim"), names(Arg), Arg)
    xArg <- as.character(Arg[c(1, m)])[2]
    yArg <- as.character(Arg[c(1, m)])[3]

    # Plot:
    if(xArg == "NULL" & yArg == "NULL") {
        yLim <- range(fullFrontier[, 2])
        xRange <- range(fullFrontier[, 1])
        xDiff <- diff(xRange)
        xLim <- c(xRange[1] - 2.5*xDiff/10, xRange[2] + xDiff/10)

        # Plot:
        if(!add){
            if(frontier == "upper" | frontier == "both") {
                plot(upperFrontier, col = col[1], xlim = xLim, ylim = yLim,
                    ann = FALSE, ...)
            } else {
                if( frontier == "both") {
                    points(fullFrontier, col = col[2],
                        xlim = xLim, ylim = yLim, ...)
                }
                if(frontier == "lower" ) {
                    plot(lowerFrontier, col = col[2],
                        xlim = xLim, ylim = yLim, ann = FALSE, ...)
                }
            }
        }
        if(frontier == "upper" | frontier == "both") {
            points(upperFrontier, col = col[1],  ...)
        }
        if(frontier == "lower" | frontier == "both") {
            points(lowerFrontier, col = col[2], ...)
        }
    } else if (xArg != "NULL" & yArg == "NULL") {
        # In this case only xlim is specified in the argument list
        yLim = range(fullFrontier[, 2])
        # Plot:
        if(!add){
            if(frontier == "upper" | frontier == "both") {
                plot(upperFrontier, col = col[1], ylim = yLim,
                    ann = FALSE, ...)
            } else {
                if( frontier == "both") {
                    points(fullFrontier, col = col[2], ylim = yLim, ...)
                }
                if(frontier == "lower" ) {
                    plot(fullFrontier, col = col[2], ylim = yLim,
                        ann = FALSE, ...)
                }
            }
        }
        if(frontier == "upper" | frontier == "both") {
            points(upperFrontier, col = col[1], ...)
        }
        if(frontier == "lower" | frontier == "both") {
            points(lowerFrontier, col = col[2], ...)
        }
    } else if(xArg == "NULL" & yArg != "NULL") {
        # In this only ylim is specified in the argument list
        xRange = range(fullFrontier[, 1])
        xDiff = diff(xRange)
        xLim = c(xRange[1] - 2.5*xDiff/10, xRange[2] + xDiff/10)
        # Plot:
        if(!add){
            if(frontier == "upper" | frontier == "both") {
                plot(upperFrontier, col = col[1], xlim = xLim,
                    ann = FALSE, ...)
            } else {
                if( frontier == "both") {
                    points(fullFrontier, col = col[2], xlim = xLim, ...)
                }
                if(frontier == "lower" ) {
                    plot(lowerFrontier, col = col[2], xlim = xLim,
                        ann = FALSE,...)
                }
            }
        }
        if(frontier == "upper" | frontier == "both") {
            points(upperFrontier, col = col[1], ...)
        }
        if(frontier == "lower" | frontier == "both") {
            points(lowerFrontier, col = col[2], ...)
        }
    } else if (xArg != "NULL" & yArg != "NULL"){
        #  If both xlim and ylim are not defined in argument list ...
        if(!add){
            if(frontier == "upper" | frontier == "both") {
                plot(fullFrontier, type = "n", ann = FALSE, ...)
                points(upperFrontier, col = col[1], ...)
            }
            if(frontier == "both") {
                points(lowerFrontier, col = col[2], ...)
            }
            if(frontier == "lower") {
                plot(lowerFrontier, col = col[2], ann = FALSE, ...)
            }
        } else{
            if(frontier == "upper" | frontier == "both") {
                points(upperFrontier, col = col[1], ...)
            }
            if(frontier == "lower" | frontier == "both") {
                points(lowerFrontier, col = col[2], ...)
            }
        }
    }

    # Add Title:
    if (title) {
        labs = attr(fullFrontier, "control")
        title(
            main = "Efficient Frontier",
            xlab = paste("Target Risk[", labs[1], "]", sep = ""),
            ylab = paste("Target Return[", labs[2], "]", sep = ""))
    }
    
    # Add Rmetrics - Do not Remove!
    mtext("Rmetrics", adj=0, side=4, cex=0.7, col="darkgrey")

    # Return Value:
    invisible(fullFrontier)
}


# ------------------------------------------------------------------------------


minvariancePoints <-
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds the minimum risk point to a MV and CVaR portfolio plot

    # Arguments:
    
    # FUNCTION:

    # Match Arguments:
    return <- match.arg(return)
    risk <- match.arg(risk)

    # Get Portfolio Slots:
    data <- getSeries(object)
    spec <- getSpec(object)
    constraints <- getConstraints(object)

    # Add Minimum Variance Point:
    mvPortfolio <- minvariancePortfolio(data, spec, constraints)
    assets <- frontierPoints(mvPortfolio, return = return, risk = risk,
        auto = auto)
    points(assets, ...)

    # Return Value:
    invisible(assets)
}


# ------------------------------------------------------------------------------


cmlPoints <-
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds the capital market line to a portfolio plot

    # Arguments:
    
    # FUNCTION:

    # Match Arguments:
    return <- match.arg(return)
    risk <- match.arg(risk)

    # Get Portfolio Statistics:
    data <- getSeries(object)
    spec <- getSpec(object)
    constraints <- getConstraints(object)

    # Add Capital Market Line Tangency Point:
    cmlPortfolio <- tangencyPortfolio(data, spec, constraints)
    assets <- frontierPoints(cmlPortfolio, return = return, risk = risk,
        auto = auto)
    points(assets, ...)

    # Return Value:
    invisible(assets)
}


# ------------------------------------------------------------------------------


cmlLines <-
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds the capital market line to a portfolio plot

    # Arguments:
    
    # FUNCTION:

    # Match Arguments:
    return <- match.arg(return)
    risk <- match.arg(risk)

    # Get Portfolio Statistics:
    data <- getSeries(object)
    spec <- getSpec(object)
    constraints <- getConstraints(object)

    # Add Capital Market Line:
    cmlPortfolio <- tangencyPortfolio(data, spec, constraints)
    riskFreeRate <- getRiskFreeRate(spec)
    slope <- ((getTargetReturn(cmlPortfolio)[, "mean"] - riskFreeRate) /
        getTargetRisk(cmlPortfolio@portfolio)[, "Cov"])
    if(slope > 0) { 
        abline(riskFreeRate, slope, ...)
    } else {
        warning("CML Line does not exist")
    }
    
    # Return Value:
    invisible(slope)
}


# ------------------------------------------------------------------------------


tangencyPoints <-
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds tangency point and line to a MV and CVaR portfolio plot

    # Arguments:
    
    # FUNCTION:

    # Match Arguments:
    return <- match.arg(return)
    risk <- match.arg(risk)

    # Get Portfolio Slots:
    data <- getSeries(object)
    spec <- getSpec(object)
    constraints <- getConstraints(object)

    # Compute Tangency Portfolio:
    tgPortfolio <- tangencyPortfolio(data, spec, constraints)

    # Add Tangency Point:
    assets <- frontierPoints(tgPortfolio, return = return, risk = risk,
        auto = auto)
    points(assets, ...)

    # Return Value:
    invisible(assets)
}


# ------------------------------------------------------------------------------


tangencyLines <-
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds tangency point and line to a MV and CVaR portfolio plot

    # Arguments:
    
    # FUNCTION:

    # Match Arguments:
    return = match.arg(return)
    risk = match.arg(risk)

    # Get Portfolio Slots:
    data <- getSeries(object)
    spec <- getSpec(object)
    constraints <- getConstraints(object)
    riskFreeRate <- getRiskFreeRate(object)

    # Compute Tangency Portfolio:
    tgPortfolio = tangencyPortfolio(data, spec, constraints)

    # Add Tangency Line:
    assets <- frontierPoints(tgPortfolio, return = return, risk = risk,
        auto = auto)
    slope <-( assets[2] - riskFreeRate ) / assets[1]
    if (slope > 0) {
        abline(riskFreeRate, slope, ...)
    } else {
        warning("Tangency point does not exist")
    }

    # Return Value:
    invisible(list(slope = slope, assets = assets))
}


# ------------------------------------------------------------------------------


equalWeightsPoints =
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds equal weights portfolio to a portfolio plot

    # Arguments:
    
    # FUNCTION:

    # Match Arguments:
    return <- match.arg(return)
    risk <- match.arg(risk)

    # Get Portfolio Statistics:
    data <- getSeries(object)
    spec <- getSpec(object)
    constraints <- getConstraints(object)
    numberOfAssets <- getNAssets(object)

    # Set Equal Weights:
    setWeights(spec) <- rep(1/numberOfAssets, times = numberOfAssets)

    # Add Equal Weights Portfolio:
    ewPortfolio <- feasiblePortfolio(data, spec, constraints)
    assets <- frontierPoints(ewPortfolio, return = return, risk = risk,
        auto = auto)
    points(assets, ...)

    # Return Value:
    invisible(assets)
}


# ------------------------------------------------------------------------------


singleAssetPoints <-
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds all single assets returns and risks to a portfolio plot

    # Arguments:
    
    # FUNCTION:

    # Add Single Assets:
    Statistics <- getStatistics(object)
    Type <- getType(object)

    # Match Arguments:
    return <- match.arg(return)
    risk <- match.arg(risk)

    # Get auto Risk:
    if (auto) {
        return = "mu"
        Type = getType(object)
        Estimator = getEstimator(object)
        if (Type == "MV") risk = "Cov"
        if (Type == "MV" & Estimator != "covEstimator") risk = "Sigma"
        if (Type == "QLPM") risk = "Sigma"
        if (Type == "CVaR") risk = "CVaR"
    }

    # Extract Return:
    if (return == "mean") {
        Return = Statistics$mean
    } else if (return == "mu") {
        Return = Statistics$mu
    }

    # Extract Risk:
    if (risk == "Cov") {
        Risk = sqrt(diag(Statistics$Cov))
    } else if (risk == "Sigma") {
        Risk = sqrt(diag(Statistics$Sigma))
    } else if (risk == "CVaR") {
        nAssets = getNAssets(object)
        Data = getSeries(object)
        alpha = getAlpha(object)
        Risk = NULL
        for (i in 1:nAssets) Risk = c(Risk, -.cvarRisk(Data[ ,i], 1, alpha))
    } else if (risk == "VaR") {
        nAssets = getNAssets(object)
        Data = getSeries(object)
        alpha = getAlpha(object)
        Risk = NULL
        for (i in 1:nAssets) Risk = c(Risk, -.varRisk(Data[ ,i], 1, alpha))
    }
    Risk = as.vector(Risk)

    # Add Points:
    assets = cbind(targetRisk = Risk, targetReturn = Return)
    attr(assets, "control") <-
        c(targetRisk = risk, targetReturn = return, auto = as.character(auto))
    points(assets, ...)

    # Return Value:
    invisible(assets)
}


# ------------------------------------------------------------------------------


twoAssetsLines <-
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds efficient long-only frontier of all portfolio pairs

    # Arguments:
    
    # Note:
    #   Only supported for "Short" and "LongOnly" Constraints!

    # FUNCTION:

    # Supported ?
    check <- rev(attr(object@constraints, "model"))[1]

    # Match Arguments:
    return <- match.arg(return)
    risk <- match.arg(risk)

    # Get Portfolio Statistics:
    data <- getSeries(object)
    Data <- getData(object)
    mu <- getMu(Data)
    Sigma <- getSigma(Data)
    spec <- getSpec(object)
    constraints <- getConstraints(object)

    # Add Frontiers for all Two-Assets Portfolios:
    N <- getNAssets(object)
    setWeights(spec) = NULL
    for (i in 1:(N-1) ) {
        for (j in (i+1):N ) {
            index = c(i, j)
            data2 = data[, index]
            Data2 = portfolioData(data2, spec)
            Data2@statistics$mu <- mu[index]
            Data2@statistics$Sigma <- Sigma[index, index]
            ans = portfolioFrontier(data = Data2, spec = spec)
            lines(frontierPoints(ans,
                return = return, risk = risk, auto = auto), ...)
        }
    }

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


sharpeRatioLines <-
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Adds Sharpe Ratio

    # Arguments:
    
    # FUNCTION:

    # Match Arguments:
    return = match.arg(return)
    risk = match.arg(risk)

    # Get Portfolio Slots:
    data <- getSeries(object)
    spec <- getSpec(object)
    constraints <- getConstraints(object)
    riskFreeRate <- getRiskFreeRate(object)
    Type <- getType(object)

    # Efficient Frontier:
    frontPoints <- frontierPoints(object, frontier = "upper",
        return = return, risk = risk, auto = auto)
    x <- frontPoints[, 1]
    y <- frontPoints[, 2] - riskFreeRate

    # Tangency Portfolio:
    tangencyPortfolio <- tangencyPortfolio(data, spec, constraints)
    # x.tg = getTargetReturn(tangencyPortfolio@portfolio)["mean"]
    x.tg = frontierPoints(tangencyPortfolio,
        return = return, risk = risk, auto = auto)[, 2]

    # Normalization to fit in EF Plot:
    norm <- x.tg / max(y/x)
    index <- 2:length(x)
    #index = index[diff(x) > 0]
    x <- x[index]
    y <- y[index]
    y.norm <- (y/x*norm)
    assets <- cbind(x, y.norm)
    lines(x, y.norm, ...)
    
    # Search for Maximum:
    index <- which.max(y.norm)
    points(x[index], y.norm[index], col = "cyan", cex = 1.5)

    # Add Tailored Labels -  2 may be a good Number ...
    x.tg <- x.tg[index]
    norm2 <- x.tg / max(y)
    Range <- range(y/x * norm)

    # Take a reasonable number of significant digits to plot, e.g. 2 ...
    nPrecision <- 3 
    Labels <- signif(Range, nPrecision)
    axis(4, at = Range, labels = c(" ", " "), cex.axis = 0.75)
    axis(4, at = mean(Range), labels = paste(Labels[1], "   ", Labels[2]),
        cex.axis = 0.75)

    # Add Axis Labels and Title:
    mtext("Sharpe Ratio", side = 4, line = 2, cex = 0.75)
    
    # Return Value:
    invisible(assets)
}


# ------------------------------------------------------------------------------


monteCarloPoints <-
function(object, mcSteps = 5000,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    auto = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Adds randomly feasible portfolios to a plot

    # Arguments:
    
    # FUNCTION:

    # Match Arguments:
    return = match.arg(return)
    risk = match.arg(risk)

    # Get Portfolio Statistics:
    Statistics = getStatistics(object)
    Type = getType(object)
    mu = Statistics$mu
    Sigma = Statistics$Sigma
    N = length(mu)

    # Get Specification:
    if (Type == "MV") {
        # Get Constraints Model:
        Model = rev(attr(object@constraints, "model"))[1]
        Model = "LongOnly"
        if (Model == "Short" | any(getConstraints(object) == "Short")) {
            # Monte Carlo Loop - Short:
            for (k in 1:mcSteps) {
                s = sign(rnorm(N, mean = rnorm(1)))
                weights = s * abs(rcauchy(N))
                weights = weights / sum(weights)
                Return = as.numeric(mu %*% weights)
                Risk = sqrt( as.numeric( t(weights) %*% Sigma %*% (weights) ) )
                points(Risk, Return, ...)
            }
        } else if (Model == "LongOnly" | any(getConstraints(object) == "LongOnly")) {
            # Monte Carlo Loop - Long Only:
            for (k in 1:mcSteps) {
                weights = abs(rcauchy(N))
                weights = weights / sum(weights)
                Return = as.numeric(mu %*% weights)
                Risk = sqrt( as.numeric( t(weights) %*% Sigma %*% (weights) ) )
                points(Risk, Return, ...)
            }
        } else {
            cat("\n\tOnly for Short and LongOnly Portfolios\n")
        }
    } else if (Type == "CVaR") {
        # Monte Carlo Loop - Long Only:
        x = getSeries(object)
        alpha = getAlpha(object)
        for (k in 1:mcSteps) {
            weights = abs(rcauchy(N))
            weights = weights / sum(weights)
            Return = as.numeric(mu %*% weights)
            Risk = .cvarRisk(x, weights, alpha)
            points(-Risk, Return, ...)
        }
    }

    # Return Value:
    invisible()
}


################################################################################


frontierPlotControl <-
function(

    # Colors:
    sharpeRatio.col   = "blue",
    minvariance.col   = "red",
    tangency.col      = "steelblue",
    cml.col           = "green",
    equalWeights.col  = "blue",
    singleAsset.col   = "topo.colors",
    twoAssets.col     = "grey",
    monteCarlo.col    = "black",

    # Point Sizes:
    minvariance.cex   = 1.25,
    tangency.cex      = 1.25,
    cml.cex           = 1.25,
    equalWeights.cex  = 1.25,
    singleAsset.cex   = 1.25,
    twoAssets.cex     = 0.01,
    monteCarlo.cex    = 0.01,
    sharpeRatio.cex   = 0.1,

    # Limits:
    xlim              = NULL,
    ylim              = NULL,

    # MC Steps:
    mcSteps           = 5000,

    # Pie Settings:
    pieR              = NULL,
    piePos            = NULL,
    pieOffset         = NULL
    )
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets frontier plot control parameters

    # Arguments:
    
    # FUNCTION:

    # Return Value:
    list(

        # Colors:
        sharpeRatio.col  = sharpeRatio.col,
        minvariance.col  = minvariance.col,
        tangency.col     = tangency.col,
        cml.col          = cml.col,
        equalWeights.col = equalWeights.col,
        singleAsset.col  = singleAsset.col,
        twoAssets.col    = twoAssets.col,
        monteCarlo.col   = monteCarlo.col,

        # Point Sizes:
        minvariance.cex  = minvariance.cex,
        tangency.cex     = tangency.cex,
        cml.cex          = cml.cex,
        equalWeights.cex = equalWeights.cex,
        singleAsset.cex  = singleAsset.cex ,
        twoAssets.cex    = twoAssets.cex,
        monteCarlo.cex   = monteCarlo.cex,
        sharpeRatio.cex  = sharpeRatio.cex,

        # Limits:
        xlim             = xlim,
        ylim             = ylim,

        # MC Steps:
        mcSteps          = 5000,

        # Pie Settings:
        pieR             = pieR,
        piePos           = piePos,
        pieOffset        = pieOffset

        )

}


################################################################################


tailoredFrontierPlot <-
function(object,
    return = c("mean", "mu"), risk = c("Cov", "Sigma", "CVaR", "VaR"),
    mText = NULL, col = NULL, xlim = NULL, ylim = NULL, 
    twoAssets = FALSE, sharpeRatio = TRUE, 
    title = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates an easy to use tailored frontier plot
    
    # Arguments:
    #   object - a portfolio object
    #   mtext - not used
    
    # FUNCTION:
    
    # 1. Plot the Frontier, add margin text, grid and ablines:
    offset <- 0.10
    risk <- match.arg(risk)
    return <- match.arg(return)
    
    # x - Range:
    if (is.null(xlim)) {
        if (risk == "Cov") {
            xmax <- max(sqrt(diag(getCov(object))))
        }
        if (risk == "Sigma") {
            xmax <- max(sqrt(diag(getSigma(object))))
        }
        if (risk == "CVaR") {
            alpha <- getAlpha(object)
            quantiles <- colQuantiles(getSeries(object), prob = alpha)
            n.max <- which.max(-quantiles)
            r <- getSeries(object)[, n.max]
            r <- r[r < quantiles[n.max]]
            xmax <- -mean(r)
        }
        if (risk == "VaR") {
            xmax <- max(-colQuantiles(getSeries(object), prob = alpha))
        }
        xlim <- c(0, xmax)
        Xlim <- c(xlim[1]-diff(xlim)*offset, xlim[2]+diff(xlim)*offset)
    } else {
        Xlim <- xlim
    }
    
    # y - Range:
    if (is.null(ylim)) {
        if (return == "mean") {
            ylim <- range(getMean(object))
        } else {
            ylim <- range(getMu(object))
        }
        Ylim <- c(ylim[1]-diff(ylim)*offset, ylim[2]+diff(ylim)*offset)
    } else {
        Ylim = ylim
    }
    
    # Frontier Plot:
    frontierPlot(object, labels = FALSE,
        return = return, risk = risk, auto = FALSE, 
        xlim = Xlim, ylim = Ylim, title = title, pch = 19, ...)
    
    # Add Grid:
    grid()
    
    # Add Center-Hair Cut:
    abline(h = 0, col = "grey")
    abline(v = 0, col = "grey")

    # 2. Add minimum risk (variance) Portfolio Point:
    data <- getData(object)
    spec <- getSpec(object)
    constraints <- getConstraints(object)
    mvPortfolio <- minvariancePortfolio(data, spec, constraints)
    minvariancePoints(object, return = return, risk = risk, auto = FALSE,
        pch = 19, col = "red")

    # 3. Add Tangency Portfolio Point and Tangency Line:
    tangencyPoints(object, return = return, risk = risk, auto = FALSE,
        pch = 19, col = "blue")
    tangencyLines(object, return = return, risk = risk, auto = FALSE,
        col = "blue")

    # 4. Add Equal Weights Portfolio:
    xy <- equalWeightsPoints(object, return = return, risk = risk, 
        auto = FALSE, pch = 15, col = "grey")
    text(xy[, 1]+diff(xlim)/20, xy[, 2]+diff(ylim)/20, "EWP",
        font = 2, cex = 0.7)

    # 5. Add all Assets Points:
    if (is.null(col)) col = rainbow(6)
    xy <- singleAssetPoints(object, return = return, risk = risk, 
        auto = FALSE, cex = 1.5,
        col = col, lwd = 2)
    text(xy[, 1]+diff(xlim)/20, xy[, 2]+diff(ylim)/20,
        rownames(xy), font = 2, cex = 0.7)

    # 6. Add optionally all Two Assets  Lines
    if (twoAssets) {
        twoAssetsLines(object, return = return, risk = risk, auto = FALSE,
            lty = 3, col = "grey")
    }

    # 6. Add Sharpe Ratio Line:
    if(sharpeRatio) {
        sharpeRatioLines(object, return = return, risk = risk, auto = FALSE,
            col = "orange", lwd = 2)
    }
    
    # Add Rmetrics - Do not Remove!
    mtext("Rmetrics", adj=0, side=4, cex=0.7, col="darkgrey")

    # Return Value:
    invisible(list(object=object, xlim=Xlim, ylim=Ylim))
}


################################################################################

