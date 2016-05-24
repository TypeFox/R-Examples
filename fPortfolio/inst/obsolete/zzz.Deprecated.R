
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                   DESCRIPTION:
#  .mvnormFit                  Fits a multivariate Normal distribution
#  .mvsnormFit                 Fits a multivariate skew-Normal distribution
#  .mvstFit                    Fits a multivariate skew-Student-t distribution
# FUNCTION:                   DESCRIPTION:
#  .assetsStats                Computes statistics of monthly assets sets 
# FUNCTION:                   DESCRIPTION:
#  .dutchPortfolioData         Example Data from Engel's Diploma Thesis
#  .usPortfolioData            Annual US Economics Portfolio Data
#  .sm132PortfolioData         Example from Scherer, Martin: Chapter 1.32
#  .worldIndexData             A data set of World Indexes
# FUNCTION:                   DESCRIPTION:
#  fixBinHistogram             Returns histogram with fixed bins
################################################################################


.mvnormFit <- 
    function(x, title=NULL, description=NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits a multivariate Normal distribution

    # Arguments:
    #   x - A multivariate time series, a data frame, or any other
    #       rectangular object of assets which can be converted into
    #       a matrix by the function as.matrix. Optional Dates are
    #       rownames, instrument names are column names.

    # Value:
    #   The function returns a list with the following entries:
    #   mu - Mean values of each asset time series
    #   Omega - Covariance matrix of assets

    # Notes:
    #   Requires function "msn.mle" ans "mst.mle" from R's GPL licensed
    #     contributed package "sn", (C) 1998-2004 A. Azzalini.
    #   The list returned by this function can serve as input for the
    #     function assetsSim().

    # FUNCTION:

    # Settings:
    assets = as.matrix(x)
    method = method[1]
    colNames = colnames(x)

    # Fit mvNormal:
    fit = list()
    mu = apply(assets, 2, mean)
    Omega = cov(assets)
    alpha = rep(0, times = length(mu))
    df = Inf

    # Add Names:
    names(mu) = colNames
    names(alpha) = colNames
    rownames(Omega) = colNames
    colnames(Omega) = colNames

    # Add Title:
    if (is.null(title))
        title = paste("Fitted Asset Data Model: ", method)

    # Add Description:
    if (is.null(description))
        description = description()

    # Return Value:
    new("fASSETS",
        call = as.call(match.call()),
        method = as.character(method),
        model = list(mu = mu, Omega = Omega, alpha = alpha, df = df),
        data = as.data.frame(x),
        fit = as.list(fit),
        title = as.character(title),
        description = as.character(description)
        )
}


# ------------------------------------------------------------------------------


.mvsnormFit <-
    function(x, title=NULL, description=NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits a multivariate skew-Normal distribution

    # Arguments:
    #   x - A multivariate time series, a data frame, or any other
    #       rectangular object of assets which can be converted into
    #       a matrix by the function as.matrix. Optional Dates are
    #       rownames, instrument names are column names.

    # Value:
    #   The function returns a list with the following entries:
    #   mu - Mean values of each asset time series
    #   Omega - Covariance matrix of assets
    #   alpha - Skewness vector

    # Notes:
    #   Requires function "msn.mle" ans "mst.mle" from R's GPL licensed
    #     contributed package "sn", (C) 1998-2004 A. Azzalini.
    #   The list returned by this function can serve as input for the
    #     function assetsSim().

    # FUNCTION:

    # Settings:
    assets = as.matrix(x)
    method = method[1]
    colNames = colnames(x)

    # Fit skew-Normal:
    fit = mvFit(assets, method = "snorm", ...)
    mu = as.vector(fit@fit$dp$beta)
    Omega = fit@fit$dp$Omega
    alpha = as.vector(fit@fit$dp$alpha)
    df = Inf
    fit = fit@fit

    # Add Names:
    names(mu) = colNames
    names(alpha) = colNames
    rownames(Omega) = colNames
    colnames(Omega) = colNames

    # Add Title:
    if (is.null(title))
        title = paste("Fitted Asset Data Model: ", method)

    # Add Description:
    if (is.null(description))
        description = description()

    # Return Value:
    new("fASSETS",
        call = as.call(match.call()),
        method = as.character(method),
        model = list(mu = mu, Omega = Omega, alpha = alpha, df = df),
        data = as.data.frame(x),
        fit = as.list(fit),
        title = as.character(title),
        description = as.character(description)
        )
}


# ------------------------------------------------------------------------------


.mvstFit <- 
    function(x, title = NULL, description=NULL, fixed.df=NA, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits a multivariate skew-Student-t distribution

    # Arguments:
    #   x - A multivariate time series, a data frame, or any other
    #       rectangular object of assets which can be converted into
    #       a matrix by the function as.matrix. Optional Dates are
    #       rownames, instrument names are column names.

    # Value:
    #   The function returns a list with the following entries:
    #   mu - Mean values of each asset time series
    #   Omega - Covariance matrix of assets
    #   alpha - Skewness vector
    #   df - Degrees of freedom, measures kurtosis

    # Notes:
    #   Requires function "msn.mle" ans "mst.mle" from R's GPL licensed
    #     contributed package "sn", (C) 1998-2004 A. Azzalini.
    #   The list returned by this function can serve as input for the
    #     function assetsSim().

    # FUNCTION:

    # Settings:
    assets = as.matrix(x)
    method = method[1]
    colNames = colnames(x)

    # Fit skew-Student:
    fit = mvFit(assets, method = "st", fixed.df = fixed.df, ...)
    mu = as.vector(fit@fit$beta)
    Omega = fit@fit$dp$Omega
    alpha = as.vector(fit@fit$dp$alpha)
    df = fit@fit$dp$df
    fit = fit@fit

    # Add Names:
    names(mu) = colNames
    names(alpha) = colNames
    rownames(Omega) = colNames
    colnames(Omega) = colNames

    # Add Title:
    if (is.null(title))
        title = paste("Fitted Asset Data Model: ", method)

    # Add Description:
    if (is.null(description))
        description = description()

    # Return Value:
    new("fASSETS",
        call = as.call(match.call()),
        method = as.character(method),
        model = list(mu = mu, Omega = Omega, alpha = alpha, df = df),
        data = as.data.frame(x),
        fit = as.list(fit),
        title = as.character(title),
        description = as.character(description)
        )
}


################################################################################


.hclustSelect <-  
    function(x, control = NULL, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description: 
    #   Hierarchical Clustering
    
    # FUNCTION:
    
    # Method:
    if (is.null(control)) 
        control = c(measure = "euclidean", method = "complete")
    measure = control[1]
    method = control[2]
    
    # hclust:
    ans = hclust(dist(t(x), method = measure), method = method, ...)
    class(ans) = c("list", "hclust")
    
    # Return Value:
    ans
}


# -----------------------------------------------------------------------------


.kmeansSelect <-  
    function(x, control = NULL, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description: 
    #   kmeans Clustering
    
    # Note:
    #   centers must be specified by the user!
    
    # FUNCTION:
    
    # Method:
    if (is.null(control)) 
        control = c(centers = 5, algorithm = "Hartigan-Wong")
    centers = as.integer(control[1])
    algorithm = control[2]
    
    # kmeans:
    ans = kmeans(x = t(x), centers = centers, algorithm = algorithm, ...)
    class(ans) = c("list", "kmeans")
    
    # Return Value:
    ans
}


################################################################################


.assetsStats <- 
    function(x)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes benchmark statistics for a data set of assets with
    #   monthly data records. 
    
    # Details:
    #   The computed statistics values are:
    #       records - number of records (length of time series)
    #       paMean - annualized (pa, per annum) Mean of Returns
    #       paAve - annualized Average of Returns
    #       paVola - annualized Volatility (standard Deviation)
    #       paSkew - Skewness of Returns
    #       paKurt - Kurtosis of Returns
    #       maxDD - maximum Drawdown
    #       TUW - Time under Water
    #       mMaxLoss - Monthly maximum Loss
    #       mVaR - Monthly 99% Value-at-Risk 
    #       mModVaR - Monthly 99% Modified Value-at-Risk 
    #       mSharpe - Monthly Sharpe Ratio
    #       mModSharpe - Monthly Modified Sharpe Ratio
    #       skPrice - Skewness/Kurtosis Price
    #   The statistics are implemented based on the formulas from
    #   "Extreme Metrics". They reflect risk measures as used in 
    #   the hedge fund software from "www.AlternativeSoft.com".

    # Arguments:
    #   x - asset data set, a matrix (or vector) where the rows
    #       are numbered by "time", and the columns belong to the
    #       individual assets. Monthly values are expected.
    
    # Value:
    #   The function returns a data frame with the values of the
    #   12 statistics for each asset.
    
    # Reference:
    #   "ExtremeMetrics Software", Help Document, Alternative Software,
    #   March 2003, 4 pages.
    
    # Example:
    
    # FUNCTION:
    
    # If x is a vector, make it a matrix:
    statistics = 14
    if (is.null(dim(x))) {
        n = 1 
        x = matrix(x, length(x)) 
        result = matrix(rep(0, times = statistics), ncol = 1) }
    else {
        n = dim(x)[2] 
        result = matrix(rep(0, times = statistics*n), ncol = n) }
    
    # Give Names to Result Matrix:  
    stat.names = c(
        "Records",      "paMean",   "paAve",    "paVola",
        "paSkew",       "paKurt",   "maxDD",    "TUW",
        "mMaxLoss",     "mVaR",     "mModVaR",  "mSharpe",
        "mModSharpe",   "skPrice")
    dimnames(result) = list(stat.names, dimnames(x)[[2]])   

    # Loop over all Assets:
    for (i in 1:n) {
        r = x[, i]
        # Number of Records:
        result[1, i] = length(r)
        # Annualized mean from monthly returns:
        result[2, i] = annualizedMean = (1 + mean(r))^12 - 1
        # Annualized mean from monthly returns:
        result[3, i] = annualizedAverage = mean(r)*sqrt(12)
        # Annualized volatility from monthly returns:
        result[4, i] = annualizedVolatility = sqrt(var(r))
        # Annualized skewness from monthly returns:
        result[5, i] = annualizedSkewness = skewness(r) 
        # Annualized Kurtosis from monthly returns:
        result[6, i] = annualizedKurtosis = kurtosis(r) 
        # Maximum Drawdown of of monthly returns:
        result[7, i] = maxDrawdown = max(cummax(cumsum(r)) - cumsum(r))
        # Time-Under-Water of monthly returns:
        result[8, i] = timeUnderWater = 
            max(diff(which (diff(cummax(cumsum(r))) != 0)))
        # Maximum Loss of monthly returns:
        result[9, i] = maxMonthlyLoss = min(r)  
        # Monthly Value at Risk:
        zc = 2.33
        result[10, i] = monthlyVaR = annualizedMean - 
            zc * annualizedVolatility   
        # Monthly Modified Value at Risk:
        p = 0.99; s = annualizedSkewness; k = annualizedKurtosis    
        zcf = zc + (zc*zc-1)*s/6 + zc*(zc*zc-3)*k/24 + zc*(2*zc*zc-5)*s*s/36
        result[11, i] = monthlyModVaR = annualizedMean - 
            zcf * annualizedVolatility  
        # Monthly Sharpe Ratio:
        result[12, i] = monthlySharpeRatio = 
            annualizedMean/annualizedVolatility 
        # Monthly Modified Sharpe Ratio:
        result[13, i] = monthlyModSharpeRatio = annualizedMean/monthlyModVaR    
        # Skewness Kurtosis Price:
        result[14, i] = skewnesskurtosisPrice = annualizedMean * 
            ( monthlyModVaR/monthlyVaR - 1) }
    
    # Result:
    ans = as.data.frame(round(result, digits = 3))    
    
    # Return Value:
    ans
} 


################################################################################


.dutchPortfolioData =
function()
{   # A function implemented by Rmetrics

    # Description:
    #   Example Portfolio Data from Engels

    # Example:
    #   engelsPortfolioData()

    # FUNCTION:

    # Mean Returns:
    mu = c(0.266, 0.274, 0.162, 0.519, 0.394, 0.231, 0.277) / 1000
    names(mu) = c(
        "Elsevier", "Fortis", "Getronics", "Heineken",
        "Philips", "RoyalDutch", "Unilever")

    # Variance-Covariance Risk:
    Sigma = c(
        0.345, 0.150, 0.183, 0.088, 0.186, 0.090, 0.095,
        0.150, 0.399, 0.204, 0.107, 0.236, 0.130, 0.127,
        0.183, 0.204, 1.754, 0.075, 0.325, 0.110, 0.091,
        0.088, 0.107, 0.075, 0.243, 0.096, 0.064, 0.086,
        0.186, 0.236, 0.325, 0.096, 0.734, 0.147, 0.114,
        0.090, 0.130, 0.110, 0.064, 0.147, 0.221, 0.093,
        0.095, 0.127, 0.091, 0.086, 0.114, 0.093, 0.219)
    Sigma = matrix(Sigma, ncol = 7)
    colnames(Sigma) = rownames(Sigma) = names(mu)

    # Return Value:
    list(mu = mu, Sigma = Sigma)
}


# ------------------------------------------------------------------------------


.usPortfolioData =
function()
{   # A function implemented by Rmetrics

    # Description:
    #   Annual US Economics Portfolio Data

    # Example:
    #   usPortfolioData()
    #   list(mu = round(mean(usPortfolioData()),5),
    #   Sigma = round(var(usPortfolioData()), 5))

    # FUNCTION:

    # Units:
    Units = c("TBills3m", "LongBonds", "SP500", "Wilshire5000",
        "NASDAQComp", "LehmanBonds", "EAFE", "Gold")

    # Time Series Object:
    tS = as.timeSeries(as.data.frame(matrix(c(
        19731231,1.075,0.942,0.852,0.815,0.698,1.023,0.851,1.677,
        19741231,1.084,1.020,0.735,0.716,0.662,1.002,0.768,1.722,
        19751231,1.061,1.056,1.371,1.385,1.318,1.123,1.354,0.760,
        19761231,1.052,1.175,1.236,1.266,1.280,1.156,1.025,0.960,
        19771231,1.055,1.002,0.926,0.974,1.093,1.030,1.181,1.200,
        19781231,1.077,0.982,1.064,1.093,1.146,1.012,1.326,1.295,
        19791231,1.109,0.978,1.184,1.256,1.307,1.023,1.048,2.212,
        19801231,1.127,0.947,1.323,1.337,1.367,1.031,1.226,1.296,
        19811231,1.156,1.003,0.949,0.963,0.990,1.073,0.977,0.688,
        19821231,1.117,1.465,1.215,1.187,1.213,1.311,0.981,1.084,
        19831231,1.092,0.985,1.224,1.235,1.217,1.080,1.237,0.872,
        19841231,1.103,1.159,1.061,1.030,0.903,1.150,1.074,0.825,
        19851231,1.080,1.366,1.316,1.326,1.333,1.213,1.562,1.006,
        19861231,1.063,1.309,1.186,1.161,1.086,1.156,1.694,1.216,
        19871231,1.061,0.925,1.052,1.023,0.959,1.023,1.246,1.244,
        19881231,1.071,1.086,1.165,1.179,1.165,1.076,1.283,0.861,
        19891231,1.087,1.212,1.316,1.292,1.204,1.142,1.105,0.977,
        19901231,1.080,1.054,0.968,0.938,0.830,1.083,0.766,0.922,
        19911231,1.057,1.193,1.304,1.342,1.594,1.161,1.121,0.958,
        19921231,1.036,1.079,1.076,1.090,1.174,1.076,0.878,0.926,
        19931231,1.031,1.217,1.100,1.113,1.162,1.110,1.326,1.146,
        19941231,1.045,0.889,1.012,0.999,0.968,0.965,1.078,0.990),
        byrow = TRUE, ncol = 9)))
    colnames(tS)<-Units

    # Return Value:
    tS
}


# ------------------------------------------------------------------------------


.sm132PortfolioData =
function()
{
    # A function implemented by Rmetrics

    # Description:
    #   Example from Scherer, Martin:  "Modern Portfolio Omtimization":
    #       Cheapter 1.32

    # FUNCTION:
    corr = matrix(data =
        c(  1, 0.4, 0.5, 0.5, 0.4, 0.1, 0.1, 0.1,
          0.4, 1.0, 0.3, 0.3, 0.1, 0.4, 0.1, 0.1,
          0.5, 0.3, 1.0, 0.7, 0.1, 0.1, 0.5, 0.1,
          0.5, 0.3, 0.7, 1.0, 0.1, 0.1, 0.1, 0.5,
          0.4, 0.1, 0.1, 0.1, 1.0, 0.0, 0.0, 0.0,
          0.1, 0.4, 0.1, 0.1, 0.0, 1.0, 0.0, 0.0,
          0.1, 0.1, 0.5, 0.1, 0.0, 0.0, 1.0, 0.2,
          0.1, 0.1, 0.1, 0.5, 0.0, 0.0, 0.2, 1.0),
          nrow = 8, ncol = 8)
    vol = diag(c(17, 21, 22, 20, 8, 8, 8, 8))
    Cov = vol %*% corr %*% vol

    # Average return
    mu = c(3, 4, 5, 6, 0.25, 0.5, 0.75, 1)

    # Return value:
    list(mu = mu, Sigma = Cov)
}


# ------------------------------------------------------------------------------


.worldIndexData =
function()
{
    # Description:
    #   A data set of World Indexs contributed by Dominik Locher

    # Units:
    Units = c("Asia", "EasternEurope", "FarEast", "LatinAmerica")

    # Time Series Object:
    x = c(
        20070327,370.04,302.41,326.56,3100.66,
        20070326,370.37,304.79,327.06,3128.91,
        20070325,369.54,302.25,326.03,3124.70,
        20070324,369.54,302.25,326.03,3124.70,
        20070323,369.54,302.25,326.03,3124.70,
        20070322,369.75,298.95,326.26,3129.17,
        20070321,365.46,292.45,322.84,3116.79,
        20070320,362.57,289.46,320.86,3034.35,
        20070319,360.93,292.24,319.81,2990.89,
        20070318,357.70,287.29,317.28,2938.57,
        20070317,357.70,287.29,317.28,2938.57,
        20070316,357.70,287.29,317.28,2938.57,
        20070315,357.74,285.52,317.04,2962.38,
        20070314,353.26,281.37,312.66,2936.81,
        20070313,362.26,285.91,320.23,2930.81,
        20070312,362.09,286.35,320.47,3014.71,
        20070311,357.45,288.41,315.81,3004.10,
        20070310,357.45,288.41,315.81,3004.10,
        20070309,357.45,288.41,315.81,3004.10,
        20070308,357.38,281.80,315.42,2964.89,
        20070307,350.68,278.35,310.37,2901.26,
        20070306,349.63,278.58,308.97,2910.81,
        20070305,342.19,273.38,302.54,2797.08,
        20070304,357.72,282.62,316.19,2880.75,
        20070303,357.72,282.62,316.19,2880.75,
        20070302,357.72,282.62,316.19,2880.75,
        20070301,359.75,280.80,317.25,2925.88,
        20070228,363.46,290.20,321.72,2957.57,
        20070227,372.72,297.04,329.05,2933.25,
        20070226,377.55,308.41,333.45,3143.55,
        20070225,378.21,304.53,334.12,3152.57,
        20070224,378.21,304.53,334.12,3152.57,
        20070223,378.21,304.53,334.12,3152.57,
        20070222,379.11,303.81,334.01,3198.17,
        20070221,378.44,300.74,332.64,3166.70,
        20070220,377.83,300.17,331.72,3157.26,
        20070219,377.94,303.03,331.21,3166.05,
        20070218,378.26,301.19,331.53,3162.13,
        20070217,378.26,301.19,331.53,3162.13,
        20070216,378.26,301.19,331.53,3162.13,
        20070215,377.28,299.89,330.64,3172.06,
        20070214,372.47,301.38,327.11,3172.37,
        20070213,368.75,295.28,323.16,3112.62,
        20070212,372.40,289.73,326.33,3049.67,
        20070211,376.56,297.99,329.20,3081.50,
        20070210,376.56,297.99,329.20,3081.50,
        20070209,376.56,297.99,329.20,3081.50,
        20070208,376.37,298.04,328.56,3111.51,
        20070207,376.14,305.12,328.39,3111.97,
        20070206,374.87,306.71,327.71,3123.29,
        20070205,372.22,304.55,324.90,3105.70,
        20070204,370.91,302.47,324.03,3096.00,
        20070203,370.91,302.47,324.03,3096.00,
        20070202,370.91,302.47,324.03,3096.00,
        20070201,366.10,302.61,319.70,3080.11,
        20070131,362.92,296.93,317.05,3041.84,
        20070130,365.45,293.86,319.34,2994.49,
        20070129,363.99,293.20,317.87,2959.63,
        20070128,365.73,295.87,319.48,3008.45,
        20070127,365.73,295.87,319.48,3008.45,
        20070126,365.73,295.87,319.48,3008.45,
        20070125,371.24,299.37,325.03,3031.37,
        20070124,372.54,298.33,326.91,3050.37,
        20070123,367.71,297.63,322.26,3005.14,
        20070122,368.07,297.03,322.01,2965.56,
        20070121,366.07,292.74,320.23,2954.21,
        20070120,366.07,292.74,320.23,2954.21,
        20070119,366.07,292.74,320.23,2954.21,
        20070118,368.51,289.85,322.62,2901.66,
        20070117,366.67,288.32,320.87,2926.80,
        20070116,367.78,292.91,322.15,2908.26,
        20070115,366.66,296.45,320.98,2933.52,
        20070114,361.66,288.98,316.46,2926.08,
        20070113,361.66,288.98,316.46,2926.08,
        20070112,361.66,288.98,316.46,2926.08,
        20070111,354.97,290.37,311.21,2902.35,
        20070110,354.90,285.22,311.93,2859.72,
        20070109,361.15,288.23,317.46,2849.87,
        20070108,362.10,304.41,318.23,2903.84,
        20070107,367.47,304.32,322.78,2880.09,
        20070106,367.47,304.32,322.78,2880.09,
        20070105,367.47,304.32,322.78,2880.09,
        20070104,370.65,307.56,325.92,2968.18,
        20070103,376.06,310.53,331.11,3002.63,
        20070102,377.21,311.52,332.33,3039.15,
        20070101,371.46,309.43,327.07,2995.67,
        20061231,371.46,309.43,327.07,2995.67,
        20061230,371.46,309.43,327.07,2995.67,
        20061229,371.46,309.43,327.07,2995.67,
        20061228,370.18,307.74,325.65,2981.90,
        20061227,368.11,304.17,323.63,2975.56,
        20061226,363.36,300.91,319.54,2926.69,
        20061225,362.36,301.54,319.41,2902.57,
        20061224,362.60,302.53,319.65,2902.57,
        20061223,362.60,302.53,319.65,2902.57,
        20061222,362.60,302.53,319.65,2902.57,
        20061221,361.54,304.50,318.98,2910.08,
        20061220,361.98,304.64,319.70,2918.35,
        20061219,356.34,300.35,313.84,2917.11,
        20061218,363.09,306.87,319.50,2936.06,
        20061217,360.37,306.83,317.06,2942.70,
        20061216,360.37,306.83,317.06,2942.70,
        20061215,360.37,306.83,317.06,2942.70,
        20061214,358.11,305.14,315.26,2938.00,
        20061213,352.99,302.33,311.23,2903.05,
        20061212,352.75,304.36,311.64,2890.34,
        20061211,356.43,305.03,314.04,2907.91,
        20061210,358.28,308.42,314.60,2895.92,
        20061209,358.28,308.42,314.60,2895.92,
        20061208,358.28,308.42,314.60,2895.92,
        20061207,363.08,308.81,318.78,2889.90,
        20061206,363.95,308.24,319.82,2891.55,
        20061205,362.05,308.20,317.71,2887.74,
        20061204,359.44,303.24,315.40,2836.95,
        20061203,360.01,300.45,316.12,2780.48,
        20061202,360.01,300.45,316.12,2780.48,
        20061201,360.01,300.45,316.12,2780.48,
        20061130,358.40,299.50,315.12,2804.62,
        20061129,354.34,296.95,311.25,2789.24,
        20061128,350.48,288.78,307.40,2726.65,
        20061127,356.66,287.56,312.97,2732.92,
        20061126,354.96,287.21,311.34,2782.22,
        20061125,354.96,287.21,311.34,2782.22,
        20061124,354.96,287.21,311.34,2782.22,
        20061123,354.65,285.92,311.04,2791.43,
        20061122,353.85,284.78,310.34,2787.78,
        20061121,349.05,284.25,305.82,2767.77,
        20061120,347.46,278.95,304.95,2740.54,
        20061119,348.12,281.08,305.55,2735.42,
        20061118,348.12,281.08,305.55,2735.42,
        20061117,348.12,281.08,305.55,2735.42,
        20061116,348.96,285.75,306.06,2761.74,
        20061115,347.24,283.87,304.79,2766.90,
        20061114,346.29,284.29,303.74,2760.15,
        20061113,343.74,283.69,301.17,2721.09,
        20061112,343.78,284.11,301.32,2733.62,
        20061111,343.78,284.11,301.32,2733.62,
        20061110,343.78,284.11,301.32,2733.62,
        20061109,343.01,283.30,300.81,2750.00,
        20061108,339.56,280.29,297.78,2750.76,
        20061107,340.64,282.54,298.59,2739.01,
        20061106,337.81,277.43,295.87,2743.48,
        20061105,338.56,275.49,296.72,2687.33,
        20061104,338.56,275.49,296.72,2687.33,
        20061103,338.56,275.49,296.72,2687.33,
        20061102,336.80,272.81,295.09,2666.32,
        20061101,333.81,277.98,292.27,2673.55,
        20061031,331.68,270.93,290.50,2663.66,
        20061030,330.43,266.78,288.86,2619.44,
        20061029,331.79,274.60,290.85,2667.38,
        20061028,331.79,274.60,290.85,2667.38,
        20061027,331.79,274.60,290.85,2667.38,
        20061026,331.41,276.15,291.20,2698.02,
        20061025,329.05,275.38,289.16,2688.38,
        20061024,328.31,272.69,288.43,2668.66,
        20061023,326.76,271.96,286.86,2654.76,
        20061022,328.10,274.18,287.91,2637.77,
        20061021,328.10,274.18,287.91,2637.77,
        20061020,328.10,274.18,287.91,2637.77,
        20061019,326.66,277.17,286.40,2651.84,
        20061018,327.51,274.63,286.91,2636.09,
        20061017,328.14,270.77,287.48,2619.25,
        20061016,329.36,271.73,288.66,2649.86,
        20061015,326.89,273.78,286.76,2625.68,
        20061014,326.89,273.78,286.76,2625.68,
        20061013,326.89,273.78,286.76,2625.68,
        20061012,322.28,267.06,282.92,2579.95,
        20061011,320.70,267.86,282.06,2558.04,
        20061010,320.94,266.72,282.39,2573.41,
        20061009,319.22,268.07,280.62,2547.12,
        20061008,323.44,262.86,284.71,2530.23,
        20061007,323.44,262.86,284.71,2530.23,
        20061006,323.44,262.86,284.71,2530.23,
        20061005,323.43,265.18,284.83,2535.34,
        20061004,320.04,259.29,282.16,2505.77,
        20061003,323.99,256.38,285.66,2449.38,
        20061002,323.89,261.75,285.40,2482.37,
        20061001,322.90,260.28,284.41,2473.06,
        20060930,322.90,260.28,284.41,2473.06,
        20060929,322.90,260.28,284.41,2473.06)
    tS = as.timeSeries(data.frame(matrix(x, byrow = TRUE, ncol = 5)))
    tS = returns(rev(tS))
    colnames(tS)<-Units

    # Return Value:
    tS
}


################################################################################


.hist <- 
    function (x, nbins) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns histogram with fixed bins
    
    # FUNCTION:
    
    # Classes:
    nclass = nbins + 1
    n = length(x)
    xname = paste(deparse(substitute(x), 500), collapse = "\n")
    
    # Breaks:
    breaks = seq(min(x), max(x), length = nclass)
    nB = length(breaks)
    h = diff(breaks)
    
    # Compute Counts:
    counts = .C("bincount", as.double(x), as.integer(n), as.double(breaks), 
        as.integer(nB), counts = integer(nB - 1), right = FALSE, 
        include = TRUE, naok = FALSE, NAOK = FALSE, DUP = FALSE, 
        PACKAGE = "base")$counts
    dens = counts/(n * h)
    mids = 0.5 * (breaks[-1] + breaks[-nB])
    
    # Histogram:
    r = structure(list(breaks = breaks, counts = counts, intensities = dens, 
        density = dens, mids = mids, xname = xname, equidist = TRUE), 
        class = "histogram")
}


################################################################################

  