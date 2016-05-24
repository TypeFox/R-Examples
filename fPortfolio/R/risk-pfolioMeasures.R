
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
#  covRisk                     Computes covariance risk as standard deviation
#  varRisk                     Computes Value at Risk
#  cvarRisk                    Computes Conditional Value at Risk
# FUNCTION:                   DESCRIPTION - DEPRECATED:
#  .covRisk                    Computes Covariance Risk
#  .varRisk                    Computes Value at Risk
#  .cvarRisk                   Computes Conditional Value at Risk
# FUNCTION:                   DESCRIPTION:
#  .cfgFit                     Fits bivariate tail dependency parameter lambda
#  .lambdaTailRisk             Fits tail lambda for multivariate data
################################################################################


covRisk <-
    function(data, weights)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Covariance Risk for assets given weights

    # Arguments:
    #   data - any univariate or multivariate object which can
    #       be transformed into a matrix
    #   weights - a numeric vector, the weights vector

    # Example:
    #   data = LPP2005.RET[, 1:6]; weights = rep(1/6, times = 6)
    #   covRisk(data, weights)

    # FUNCTION:

    # Data:
    if (inherits(data, "timeSeries"))
        data <- getDataPart(data)
    Data <- as.matrix(data)
    nAssets = dim(Data)[2]

    # Covariance Matrix:
    Sigma = cov(Data)

    # Risk:
    weights = as.vector(weights)
    Std = sqrt( as.numeric( weights %*% Sigma %*% weights ) )
    names(Std) = "Cov"

    # Return Value:
    Std

}


# ------------------------------------------------------------------------------


varRisk <-
    function(data, weights, alpha = 0.05)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes VaR for assets given weights and alpha

    # Arguments:
    #   data - any univariate or multivariate object which can
    #       be transformed into a matrix
    #   weights - a numeric vector, the weights vector
    #   alpha - a numeric value, the quantile

    # Example:
    #   data = LPP2005.RET[, 1:6]; weights = rep(1/6, times = 6)
    #   varRisk(data, weights)

    # FUNCTION:
    if (inherits(data, "timeSeries"))
        data <- getDataPart(data)

    # VaR:
    weights = as.vector(weights)
    X = as.matrix(data) %*% weights
    VaR = quantile(X, alpha, type = 1)
    names(VaR) <- paste("VaR.", alpha*100, "%", sep = "")

    # Return Value:
    VaR
}


# ------------------------------------------------------------------------------


cvarRisk <-
    function(data, weights, alpha = 0.05)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes CVaR for assets given weights and alpha

    # Arguments:
    #   data - any univariate or multivariate object which can
    #       be transformed into a matrix
    #   weights - a numeric vector, the weights vector
    #   alpha - a numeric value, the quantile

    # Example:
    #   data = LPP2005.RET[, 1:6]; weights = rep(1/6, times = 6)
    #   cvarRisk(data, weights)

    # FUNCTION:
    if (inherits(data, "timeSeries"))
        data <- getDataPart(data)

    # CVaR:
    weights = as.vector(weights)
    X = as.matrix(data) %*% weights
    VaR = quantile(X, alpha, type = 1)
    CVaR = c(CVaR = VaR - 0.5 * mean(((VaR-X) + abs(VaR-X))) / alpha)
    names(CVaR) <- paste("CVaR.", alpha*100, "%", sep = "")

    # Return Value:
    CVaR
}


################################################################################
# OLD FUNCTIONS:
#   check where they are still used


.covRisk <-
    function(data, weights)
{
    # A function implemented by Rmetrics

    # Description:
    #   Computes Covariance Risk for assets given weights and alpha

    # FUNCTION:
    if (inherits(data, "timeSeries"))
        data <- getDataPart(data)

    # Data:
    Data = as.matrix(data)
    nAssets = dim(Data)[2]

    # Mean Vector and Covariance:
    mu = colMeans(Data)
    Sigma = cov(Data)

    # Return and Risk:
    return = as.numeric( weights %*% mu )
    risk = sqrt( as.numeric( weights %*% Sigma %*% weights ) )

    # Return Value:
    list(risk = risk, return = return)
}


# ------------------------------------------------------------------------------


.varRisk <-
    function(x, weights, alpha = 0.05)
{
    # A function implemented by Rmetrics

    # Description:
    #   Computes VaR for assets given weights and alpha

    # Arguments:
    #   x - any univariate or multivariate object which can
    #       be transformed into a matrix
    #   weights - a numeric vector, the weights vector
    #   alpha - a numeric value, the quantile

    # FUNCTION:
    if (inherits(x, "timeSeries"))
        x <- getDataPart(x)

    # VaR:
    X = as.matrix(x) %*% weights
    VaR = quantile(X, alpha, type = 1)
    names(VaR) <- paste("VaR.", alpha*100, "%", sep = "")

    # Return Value:
    VaR
}


# ------------------------------------------------------------------------------


.cvarRisk <-
    function(x, weights, alpha = 0.05)
{
    # A function implemented by Rmetrics

    # Description:
    #   Computes CVaR for assets given weights and alpha

    # Arguments:
    #   x - any univariate or multivariate object which can
    #       be transformed into a matrix
    #   weights - a numeric vector, the weights vector
    #   alpha - a numeric value, the quantile

    # FUNCTION:
    if (inherits(x, "timeSeries"))
        x <- getDataPart(x)

    # CVaR:
    X = as.matrix(x) %*% weights
    VaR = quantile(X, alpha, type = 1)
    CVaR = c(CVaR = VaR - 0.5 * mean(((VaR-X) + abs(VaR-X))) / alpha)
    names(CVaR) <- paste("CVaR.", alpha*100, "%", sep = "")

    # Return Value:
    CVaR
}


################################################################################


.cfgFit <-
    function(x, y, tail = c("upper", "lower"))
{
    # Description:
    #   Fits bivariate tail dependency parameter lambda

    # Arguments:
    #   data - multivariate time series object of class S4 or a numeric
    #       matrix
    #   tail - which tail should be considered?

    # FUNCTION:

    # Match Arguments:
    tail = match.arg(tail)

    # If Lower Tail:
    if(tail == "lower") {
        x = 1-x
        y = 1-y
    }

    # Fit lambda:
    lambda = NULL
    n = length(x)
    for(i in 1:n){
        lambda = c(lambda,
            log(sqrt(log(1/x[i])*log(1/y[i]))/log(1/max(x[i],y[i])^2)))
    }
    ans = (2-2*exp(sum(lambda/n)))
    attr(ans, "control") <- c(tail = tail)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.lambdaTailRisk <-
    function(data, tail = c("upper", "lower"), margins = "norm", ...)
{
    # Description:
    #   Fits tail dependency parameter lambda for multivariate data

    # Arguments:
    #   data - multivariate time series object of class S4 or a numeric
    #       matrix
    #   tail - which tail should be considered?

    # Example:
    #   r = rarchmCopula(60, alpha = 2, type = "4")
    #   .cfgFit(r[, 1], r[, 2])
    #   x = cbind(qnorm(r[, 1]), qnorm(r[, 2]))
    #   .lambdaTailRisk(x)

    # FUNCTION:

    # Check Data:
    if(is.timeSeries(data)) data = series(data)
    n = ncol(data)

    # Normal Margins - Create Bivariate Copulae:
    x = data
    for (i in 1:n) {
        y = as.vector(data[, i])
        x[, i] = pnorm(y, mean(y), sd(y))
    }

    # Match Arguments:
    tail = match.arg(tail)

    # Compute Tail Risks:
    riskMatrix = diag(n)
    # Compute lambda:
    for ( i in 1:(n-1) ) {
        for ( j in (i+1):n ) {
            riskMatrix[i, j] = riskMatrix[j, i] =
                .cfgFit(x[, i], x[, j], tail = tail)
        }
    }
    attr(riskMatrix, "control") <- c(tail = tail)

    # Return Value:
    riskMatrix
}


################################################################################

