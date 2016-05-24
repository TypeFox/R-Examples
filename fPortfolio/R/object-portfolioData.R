
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
#  portfolioData            Returns an object of class fPFOLIODATA
################################################################################


portfolioData <-
    function(data, spec=portfolioSpec())
{
    # A function implemented by Rmetrics

    # Description:
    #   Creates portfolio data list

    # Arguments:
    #   data - a multivariate 'timeSeries' object
    #   spec -  a portfolio specification structure, from which
    #       the mean and covariance type of estimator will be extracted
    
    # Details:
    #   The first argument can be either:
    #   1) an object of class "fPFOLIODATA"
    #   2) an object of class "timeSeries"
    #   3) a "list: with portfolio mean and covariance

    # FUNCTION:

    # Data, if we have already an object of class "fPFOLIODATA":
    if (is(data, "fPFOLIODATA")) return(data)

    # Data, if we have a "timeSeries" or a "list":
    if (class(data) == "timeSeries") {
        series = data = sort(data)
        assetsNames = colnames(data)
    } else if (class(data) == "list") {
        series = rep(NA, times = length(data[[1]]))
        assetsNames = names(series) = names(data[[1]])
    }
    nAssets = length(assetsNames)
    names = assetsNames
    if(is.null(names)) names = paste("A", 1:nAssets, sep = "")
    .data = list(
        series = series,
        nAssets = nAssets,
        names = assetsNames)
        
    # Statistics:
    if (class(data) == "timeSeries") {
        estimator = getEstimator(spec)
        estimatorFun = match.fun(estimator)
        muSigma = estimatorFun(data, spec)
        Cov = cov(data)
        rownames(Cov) <- colnames(Cov) <- names
        .statistics = list(
            mean = colMeans(data),
            Cov = Cov,
            estimator = estimator,
            mu = muSigma$mu,
            Sigma = muSigma$Sigma)
    } else if (class(data) == "list") {
        .statistics = list(
            mean = data[[1]],
            Cov = data[[2]],
            estimator = NA,
            mu = data[[1]],
            Sigma = data[[2]])
    }

    # Tail Risk:
    .tailRisk = spec@model$tailRisk

    # Return Value:
    new("fPFOLIODATA",
        data = .data,
        statistics = .statistics,
        tailRisk = .tailRisk)
}


################################################################################

