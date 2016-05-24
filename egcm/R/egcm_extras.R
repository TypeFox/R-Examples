# egcm_extras.R
# Copyright (C) 2014 by Matthew Clegg

# A few extra helper routines for the egcm module.

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
#  http://www.r-project.org/Licenses/

yegcm <- function(ticker1, ticker2, 
    start=as.numeric(format(Sys.Date()-365, "%Y%m%d")),   # Starting date 
    end=as.numeric(format(Sys.Date(), "%Y%m%d")),         # Ending date 
    ...) { # Additional parameters to be passed to egcm
    # Fetches the price series of two securities from Yahoo
    # and constructs a cointegration model from them.

#    require(TTR)
    p1 <- getYahooData(ticker1, start, end)
    p2 <- getYahooData(ticker2, start, end)
    prices <- cbind(p1$Close, p2$Close)
    colnames(prices) <- c(ticker1, ticker2)
	e <- egcm(prices, ...)
	e
}

allpairs.egcm <- function (tickers, # A list of ticker symbols of securites
    startdate = as.numeric(format(Sys.Date()-365, "%Y%m%d")), # Starting date
    enddate = as.numeric(format(Sys.Date(), "%Y%m%d")),       # Ending date
    ...  # Additional parameters for egcm
) {
    # Constructs Engle-Granger cointegration models for all pairs
    # of tickers, using data from startdate to enddate.  The data
    # is fetched from Yahoo! using TTR.
    #
    # Copyright (C) 2014 Matthew Clegg
    #
    # Example:
    #   allpairs.egcm(c("SPY","DIA","QQQ","VOO"), 20130101, 20131231)
    #
    # Input parameters:
    #   tickers:    A list of ticker symbols
    #   startdate:  Starting date for the data.  Given as YYYYMMDD
    #   enddate:    Ending date for the data.  Given as YYYYMMDD
    #   ...:        Additional parameters to be passed to egcm function
    #
    # Returns:
    #   A data.frame containing the results of performing pairwise
    #   cointegration tests.  The columns of the data.frame are as
    #   follows:
    #
    #	series1:  Name of the first ticker in this cointegration test
    #	series2:  Name of the second ticker in this cointegration test
    #   log:      Boolean which if TRUE indicates that the cointegration
    #             test is performed on the logs of the series
    #   i1test:   Name of the test used for checking that the series are
    #             integrated.
    #   urtest:   Name of the test used for checking for a unit root in
    #             the residual series
    #   alpha:    Constant term of the linear relation between the series
    #   alpha.se: Standard error of alpha
    #   beta:     Linear term of the linear relation between the series
    #   beta.se:  Standard error of beta
    #   rho:      Coefficient of mean reversion
    #   rho.se:   Standard error of rho
    #   s1.i1.stat: Statistic computed for integration test of first series
    #   s1.i1.p:   p-value for integration test of first series
    #   s2.i1.stat: Statistic computed for integration test of second series
    #   s2.i1.p:   p-value for integration test of second series
    #   r.stat:    Statistic computed for cointegration test (e.g. whether
    #              the residual series contains a unit root)
    #   r.p:       p-value associated with r.stat
    #   eps.ljungbox.stat:  Ljung-Box statistic computed on the innovations
    #              of the series
    #   eps.ljungbox.p:  p-value associated with the Ljung-Box statistic
    #   s1.dsd:    Standard deviation of the first differences of the first series
    #   s2.dsd:    Standard deviation of the first differences of the second series
    #   residuals.sd: Standard deviation of the residual series
    #   eps.sd:    Standard deviation of the innovations
    #   is.cointegrated:  TRUE if the pair is cointegrated at the 5% confidence level

#    require(TTR)

    if (missing(tickers) || (length(tickers) < 2)) {
        stop("tickers must contain at least two ticker symbols")
    }
    
    if (!is.null(dim(tickers)) && length(dim(tickers)) > 1) {
        data <- as.data.frame(tickers)
        tickers <- colnames(data)
    } else {
        data <- do.call("cbind", lapply(tickers, 
            function(t) try(getYahooData(t, startdate, enddate)$Close, silent=TRUE)))
        colnames(data) <- tickers
    }
    
    if (any(is.na(data))) {
        missing_cols <- apply(data, 2, function(x) any(is.na(x)))
        stop ("There were missing observations in the data for ", tickers[missing_cols])
    }
    
    idf <- expand.grid(1:length(tickers), 1:length(tickers))
    idf <- idf[idf[,1] < idf[,2],]
    egcm.df <- do.call("rbind", lapply(1:nrow(idf), 
        function(k) as.data.frame(egcm(data[,c(idf[k,1], idf[k,2])], ...))))
    egcm.df$is.cointegrated <- (egcm.df$s1.i1.p > 0.05) & (egcm.df$s2.i1.p > 0.05) &
        (egcm.df$r.p <= 0.05)
        
    egcm.df
}
