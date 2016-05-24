################################################################################
#' S&P 500: Standard and Poor's 500 stock index, 2007--2010
#'
#' Contains the returns of the S&P 500 stock index for the years 2007--2010.
#' The returns were computed as \code{(Adjusted.Close-Open)/Open}.
#'
#' The data was downloaded from the Yahoo! Finance Website.
#'
#' @format A univariate time series with 1008 observations; a \code{zoo} object
#'
#' @name data-sp500
#' @aliases sp500
#' @docType data
#'
#' @references Yahoo! Finance Website
#' \url{http://finance.yahoo.com/q/hp?s=^GSPC+Historical+Prices}
#' @keywords data
#'
#' @examples
#' plot(sp500)
################################################################################
NULL

################################################################################
#' Beveridge's Wheat Price Index (detrended and demeaned), 1500--1869
#'
#' Contains a detrended and demeaned version of the well-known Beveridge Wheat
#' Price Index which gives annual price data from 1500 to 1869, averaged over
#' many locations in western and central Europe [cf. Beveridge (1921)].
#' The index series \eqn{x_t} was detrended as proposed by Granger (1964), p. 21, by
#' letting
#' \deqn{y_t := \frac{x_t}{\sum_{j=-15}^{15} x_{t+j}},}
#' where \eqn{x_t := x_1, t < 1} and \eqn{x_t := x_n, t > n}.
#' The time series in the data set is also demeaned by letting
#' \deqn{z_t := y_t - n^{-1} \sum_{t=1}^n y_t.}
#'
#' The index data cited in Beveridge's paper was taken from \code{bev} in the
#' \code{tseries} package.
#'
#' @format A univariate time series \eqn{(z_t)} with 370 observations; a \code{ts} object.
#'
#' @name data-wheatprices
#' @aliases wheatprices
#' @docType data
#'
#' @references
#' Beveridge, W. H. (1921). Weather and Harvest Cycles. \emph{The Economic Journal},
#' 31(124):429--452.
#'
#' Granger, C. W. J. (1964). \emph{Spectral Analysis of Economic Time Series}.
#' Princeton University Press, Princeton, NJ.
#'
#' @keywords data
#'
#' @examples
#' plot(wheatprices)
################################################################################

## This script was used to yield the data:
#require(tseries)
#data(bev)
#res <- bev[1:length(bev)]
#
## Detrend
#temp <- c(rep(res[1],15), res, rep(res[370],15))
#detrendedwp <- rep(0,370)
#for (i in 1:370) {
#  detrendedwp[i] <- res[i] / sum(temp[i+15+(-15):(15)])
#}
## Demean
#wheatprices <- ts(detrendedwp-mean(detrendedwp), start=1500, end=1869)

NULL

