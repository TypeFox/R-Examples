
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
# FUNCTION:                   DESCRIPTION:
#  pfolioVaR                   Computes VaR for a portfolio of assets
#  pfolioCVaR                  Computes CVaR for a portfoluio of assets
#  pfolioCVaRplus              Computes CVaR-Plus for a portfolio of assets
#  lambdaCVaR                  Computes CVaR's atomic split value lambda
# FUNCTION:                   DESCRIPTION:
#  pfolioMaxLoss               Computes maximum loss for a portfolio 
#  pfolioReturn                Computes return series for a portfolio
#  pfolioTargetReturn          Computes target return for a portfolio
#  pfolioTargetRisk            Computes target risk for a portfolio
#  pfolioHist                  Plots a histogram of portfolio returns
################################################################################


pfolioVaR <-
function(x, weights = NULL, alpha = 0.05)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute Value-at-Risk for a portfolio of assets

    # Arguments:
    #   x - a time series, data.frame or any other rectangular object
    #       of assets which can be written as a matrix object
    #   weights - a numeric vector of weights
    #   alpha - the confidence level

    # FUNCTION:

    # Transform:
    x <- as.matrix(x)

    # Compute Portfolio VaR:
    if (is.null(weights))
        weights <- rep(1/dim(x)[[2]], dim(x)[[2]])
    n <- dim(x)[1]
    x <- apply(t(t(x) * weights), 1, sum)
    n.alpha <- max(floor(n * alpha))
    ans <- as.vector(sort(x)[n.alpha])
    names(ans) <- "VaR"

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


pfolioCVaRplus <-
function(x, weights = NULL, alpha = 0.05)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute Value-at-Risk Plus for a portfolio of assets

    # Arguments:
    #   x - a time series, data.frame or any other rectangular object
    #       of assets which can be written as a matrix object
    #   weights - a numeric vector of weights
    #   alpha - the confidence level

    # FUNCTION:

    # Transform:
    x = as.matrix(x)

    # Compute Portfolio CVaRplus:
    if (is.null(weights)) {
        weights = rep(1/dim(x)[[2]], dim(x)[[2]])
    }
    n = dim(x)[1]
    x = apply(t(t(x) * weights), 1, sum)
    n.alpha = max(1, floor(n * alpha)-1)
    ans = as.vector(mean(sort(x)[1:n.alpha]))
    names(ans) = "CVaRplus"

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


pfolioCVaR <-
function(x, weights = NULL, alpha = 0.05)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute Conditional Value-at-risk for a portfolio of assets

    # Arguments:
    #   x - a time series, data.frame or any other rectangular object
    #       of assets which can be written as a matrix object
    #   weights - a numeric vector of weights
    #   alpha - the confidence level
    #   lambda - split value

    # FUNCTION:

    # Transform:
    data = as.matrix(x)

    # Input Data:
    if (is.null(weights)) {
        weights = rep(1/dim(data)[[2]], dim(data)[[2]])
    }
    n = dim(data)[1]
    Rp = apply(t(t(data)*weights), 1, sum)

    # Sort the Portfolio returns Y
    sorted = sort(Rp)

    # Compute Portfolio VaR:
    n.alpha = floor(n*alpha)
    VaR = sorted[n.alpha]

    # Compute Portfolio CVaRplus:
    n.alpha = max(1, floor(n*alpha)-1)
    CVaRplus = mean(sorted[1:n.alpha])

    # Compute Portfolio CVaR:
    lambda = 1 - floor(n*alpha)/(n*alpha)
    ans = as.vector(lambda*VaR + (1-lambda)*CVaRplus)
    names(ans) = "CVaR"
    attr(ans, "control") = c(CVaRplus = CVaRplus, lambda = lambda)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


lambdaCVaR <-
function(n, alpha = 0.05)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes CVaR's atomic split value lambda

    # Arguments:
    #   n - the number of oberservations
    #   alpha - the confidence interval

    # FUNCTION:

    # Compute CVaR lambda:
    lambda = 1 - floor(alpha * n) / (alpha * n)
    names(lambda) = "lambda"

    # Return Value:
    lambda
}


################################################################################


pfolioMaxLoss <-
function(x, weights = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes maximum loss for a portfolio of assets

    # Arguments:
    #   x - a timeSeries, data.frame or any other rectangular object
    #       of assets which can be written as a matrix object
    #   weights - the vector of weights
    #   alpha - the confidence level

    # FUNCTION:

    # Transform:
    x = as.matrix(x)

    # Compute MaxLoss [MinReturn]:
    if (is.null(weights)) {
        weights = rep(1/dim(x)[[2]], dim(x)[[2]])
    }
    x = apply(t(t(x)*weights), 1, sum)
    ans = min(x)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


pfolioReturn <- 
  function(x, weights=NULL, geometric=FALSE)
{   
  # A function implemented by Diethelm Wuertz
  
  # Description:
  #    Returns portfolio returns
  
  # Arguments:
  #    x - a 'timeSeries' object
  
  # Details:
  #    A fast(er) reimplementation
  
  # FUNCTION:
  
  # Compute Portfolio Returns:
  weights <- as.vector(weights)
  if(geometric) {
    X <- t ( colCumprods(1+x) - 1 )
    X <- rbind(diff(  t ( X * weights ) ))
    Return <- x[, 1]
    series(Return[+1, ]) <- x[1, ] %*% weights
    series(Return[-1, ]) <- rowSums(X)
  } else {
    Return <- x[, 1]
    series(Return) <- x %*% weights
  }
  colnames(Return) <- "pfolioRet"
  
  # Return Value:
  Return
}


# -----------------------------------------------------------------------------



# ------------------------------------------------------------------------------


pfolioTargetReturn <-
function(x, weights = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes return value of a portfolio

    # Arguments:
    #   x - a timeSeries, data.frame or any other rectangular object
    #       of assets which can be written as a matrix object
    #   weights - the vector of weights

    # FUNCTION:

    # Transform:
    x = as.matrix(x)

    # Compute Portfolio Returns:
    ans = mean(pfolioReturn(x = x, weights = weights))

    # Return Value:
    names(ans) = "TargetReturn"
    ans
}


# ------------------------------------------------------------------------------


pfolioTargetRisk <-
function(x, weights = NULL)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes risk from covariance matrix of a portfolio

    # Arguments:
    #   x - a timeSeries, data.frame or any other rectangular object
    #       of assets which can be written as a matrix object
    #   weights - the vector of weights

    # FUNCTION:

    # Transform:
    x = as.matrix(x)

    # Compute Portfolio Returns:
    if (is.null(weights))
        weights = rep(1/dim(x)[[2]], dim(x)[[2]])
    ans = as.vector(sqrt(weights %*% cov(x) %*% weights))

    # Return Value:
    names(ans) = "TargetRisk"
    ans
}


# ------------------------------------------------------------------------------


pfolioHist <-
function(x, weights = NULL, alpha = 0.05, range = NULL, details = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots a histogram of the returns of a portfolio

    # Arguments:
    #   x - a timeSeries, data.frame or any other rectangular object
    #       of assets which can be written as a matrix object
    #   weights - the vector of weights

    # FUNCTION:

    # Transform:
    x = as.matrix(x)

    # Suppress Warnings:
    opt = options()
    options(warn = -1)

    # Plot Portfolio Returns:
    Returns = pfolioReturn(x = x, weights = weights)
    if (is.null(range)) {
        lim = 1.05 * pfolioMaxLoss(x = x, weights = weights)[[1]]
        xlim = c(lim, -lim)
    } else {
        xlim = range
    }
    Histogram = hist(Returns, xlim = xlim, xlab = "Portfolio Return %",
        probability = TRUE, col = "steelblue4", border = "white", ...)
    r = seq(xlim[1], xlim[2], length = 201)
    lines(r, dnorm(r, mean = mean(Returns), sd = sd(Returns)), ...)

    points(Returns, rep(0, length(Returns)), pch = 20,
        col = "orange", cex = 1.25)

    # Add VaR, CVaRplus and MaxLoss:
    V1 = pfolioVaR(x = x, weights = weights, alpha = alpha)[[1]]
    abline(v = V1, col = "blue", ...)
    V2 = pfolioCVaRplus(x = x, weights = weights, alpha = alpha)[[1]]
    abline(v = V2, col = "red", ...)
    V3 = pfolioMaxLoss(x = x, weights = weights)[[1]]
    abline(v = V3, col = "green", ...)
    V4 = as.vector(mean(Returns))[[1]]
    V5 = as.vector(sd(Returns))[[1]]

    yt = max(density(Returns)$y)

    text(V1, yt, as.character(round(V1, 2)), cex = 0.75, col = "orange")
    text(V2, yt, as.character(round(V2, 2)), cex = 0.75, col = "orange")
    text(V3, yt, as.character(round(V3, 2)), cex = 0.75, col = "orange")
    text(V4, yt, as.character(round(V4, 2)), cex = 0.75, col = "orange")

    yt = 0.95 * yt
    text(V1, yt, "VaR", cex = 0.75, col = "orange")
    text(V2, yt, "CVaR+", cex = 0.75, col = "orange")
    text(V3, yt, "maxLoss", cex = 0.75, col = "orange")
    text(V4, yt, "Mean", cex = 0.75, col = "orange")

    # Result:
    options(opt)
    ans = list(VaR = V1, VaRplus = V2, maxLoss = V3, mean = V4, sd = V5)
    if (details) {
        cat("\nVaR:      ", V1)
        cat("\nVaRplus:  ", V2)
        cat("\nmax Loss: ", V3)
        cat("\nMean:     ", V4)
        cat("\nStDev:    ", V5)
        cat("\n")
    }

    # Return Value:
    invisible(ans)
}


################################################################################

