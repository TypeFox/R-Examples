#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# FUNCTION:                 DESCRIPTION:
#  endOfPeriodSeries         Returns series back to a given period
#  endOfPeriodStats          Returns statistics back to a given period
#  endOfPeriodBenchmarks     Returns benchmarks back to a given period
################################################################################


endOfPeriodSeries <- 
    function(x, nYearsBack = c("1y", "2y", "3y", "5y", "10y", "YTD"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns series back to a given period

    # Arguments:
    #   x - a monthly 'timeSeries' object of returns
    #   nYearsBack - a period string. How long back should the series
    #       be extracted? Options include values from 1 year to 10 years.
    #       and year to date: "1y", "2y", "3y", "5y", "10y", "YTD". 

    # FUNCTION:
    
    # Check:
    stopifnot(is.timeSeries(x))
      
    # Match Arguments:
    nYearsBack <- match.arg(nYearsBack)
    
    # Settings:
    if (nYearsBack == "YTD") monthsBack = atoms(end(x))$m else 
    if (nYearsBack == "1y") monthsBack = 12 else 
    if (nYearsBack == "2y") monthsBack = 24 else 
    if (nYearsBack == "3y") monthsBack = 36 else 
    if (nYearsBack == "5y") monthsBack = 60 else 
    if (nYearsBack == "10y") monthsBack = 120
    stopifnot( nrow(x) >= monthsBack )
      
    # ReturnValue:
    rev(rev(x)[1:monthsBack, ])
}


# ------------------------------------------------------------------------------


endOfPeriodStats <- 
    function(x, nYearsBack = c("1y", "2y", "3y", "5y", "10y", "YTD"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns series statistics back to a given period

    # Arguments:
    #   x - a monthly 'timeSeries' object of returns
    #   nYearsBack - a period string. How long back should the series
    #       be extracted? Options include values from 1 year to 10 years.
    #       and year to date: "1y", "2y", "3y", "5y", "10y", "YTD".

    # FUNCTION:
    
    # Check:
    stopifnot(is.timeSeries(x))

    # Match Arguments:
    nYearsBack <- match.arg(nYearsBack)

    # Series:
    Series <- endOfPeriodSeries(x, nYearsBack = nYearsBack)

    # Internal Function:
    .cl.vals <- function(x, ci) {
        x = x[!is.na(x)]
        n = length(x)
        if (n <= 1) return(c(NA, NA))
        se.mean = sqrt(var(x)/n)
        t.val = qt((1 - ci)/2, n - 1)
        mn = mean(x)
        lcl = mn + se.mean * t.val
        ucl = mn - se.mean * t.val
        c(lcl, ucl)
    }

    # Statistics:
    for (i in 1:ncol(Series))
    {
        # Basic Statistics:
        X = as.vector(Series[, i])
        X.length = length(X)
        X = X[!is.na(X)]
        X.na = X.length - length(X)
        ci = 0.95
        z = c(X.length, X.na, min(X), max(X),
            as.numeric(quantile(X, prob = 0.25, na.rm = TRUE)),
            as.numeric(quantile(X, prob = 0.75, na.rm = TRUE)),
            mean(X), median(X), sum(x), sqrt(var(X)/length(X)),
            .cl.vals(X, ci)[1], .cl.vals(X, ci)[2],
            var(X), sqrt(var(X)), skewness(X), kurtosis(X))
        znames = c("nobs", "NAs", "Minimum", "Maximum", "1. Quartile",
            "3. Quartile", "Mean", "Median", "Sum", "SE Mean",
            "LCL Mean", "UCL Mean", "Variance", "Stdev", "Skewness",
            "Kurtosis")
        stats1 <- matrix(z, ncol = 1)
        row.names(stats1) <- znames

        # Monthly Return Statistics:
        xData <- as.vector(x)
        noNegativePeriods <- length(xData[xData < 0 ])
        noPositivePeriods <- length(xData[xData > 0 ])
        stats1 = rbind(stats1,
            worstPeriod = min(xData),
            negativeValues = noNegativePeriods,
            positiveValues = noPositivePeriods)

        MaximumDrawdown = NA
        TimeUnderWater = NA
        AnnualizedVolatility = NA
        SharpeRatio = NA
        InformationRatio = NA
        ValueAtRisk = NA
        ExpectedShortfall = NA

        # Bind:
        if (i > 1) {
            stats <- cbind.data.frame(stats, stats1)
        } else {
            stats <- stats1
        }
    }
    colnames(stats) <- colnames(x) 
    

    # Return Value:
    stats
}


# ------------------------------------------------------------------------------


endOfPeriodBenchmarks <- 
    function(x, benchmark = ncol(x),
    nYearsBack = c("1y", "2y", "3y", "5y", "10y", "YTD"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns benchmarks back to a given period

    # Arguments:
    #   x - a monthly 'timeSeries' object of financial returns
    #   nYearsBack - a period string. How long back should the series
    #       be extracted? Options include values from 1 year to 10 years.
    #       and year to date: "1y", "2y", "3y", "5y", "10y", "YTD".

    # FUNCTION:
    
    # Checks:
    stopifnot(is.timeSeries(x))

    # Match Arguments:
    nYearsBack <- match.arg(nYearsBack)

    # Series:
    Series <- endOfPeriodSeries(x[, -benchmark], nYearsBack = nYearsBack)
    y <- Benchmark <- endOfPeriodSeries(x[, benchmark], nYearsBack = nYearsBack)

    stats <- NULL
    for (i in 1:ncol(Series))
    {
        # Gdet Series:
        x <- Series[, i]

        # Compute Statistics:
        stats1 <- c(
            TrackingError = NA,
            Alpha = NA,
            Beta = NA,
            CorrelationToBenchmark = NA)

        # Bind Results:
        stats <- rbind(stats, stats1)
    }

    # Return Value:
    invisible()
}


################################################################################

