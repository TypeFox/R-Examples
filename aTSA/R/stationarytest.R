#' Stationary Test for Univariate Time Series
#' @description Performs stationary test for a univariate time series.
#' @param x a numeric vector or univariate time series.
#' @param method a character indicating which test to use. The default is 
#' \code{"adf"} by Augmented Dickey-Fuller test. 
#' @param nlag the lag order to calculate the test statistic, only valid for 
#' \code{method = "adf"}. See \code{\link{adf.test}} for more details.
#' @param type the test type, only valid for \code{method = "pp"}. 
#' See \code{\link{pp.test}} for more details.
#' @param lag.short a logical value, only valid for \code{method = "pp"} or \code{"kpss"}. 
#' See \code{\link{pp.test}} and \code{\link{kpss.test}} for more details.
#' @param output a logical value indicating to print the results in R console. 
#' The default is \code{TRUE}.
#' @details This function combines the existing functions \code{\link{adf.test}}, 
#' \code{\link{pp.test}} and 
#' \code{\link{kpss.test}} for testing the stationarity of a univariate time series \code{x}.
#' @return The results are the same as one of the \code{\link{adf.test}}, \code{\link{pp.test}},
#' \code{\link{kpss.test}}, depending on which test are used.
#' 
#' @note Missing values are removed.
#' @author Debin Qiu
#' @examples x <- arima.sim(list(order = c(1,0,0),ar = 0.2),n = 100)
#' stationary.test(x)  # same as adf.test(x)
#' stationary.test(x, method = "pp") # same as pp.test(x)
#' stationary.test(x, method = "kpss") # same as kpss.test(x)
#' 
#' @export 
stationary.test <- function(x,method = c("adf","pp","kpss"),nlag = NULL,
                            type = c("Z_rho","Z_tau"),
                            lag.short = TRUE, output = TRUE)
{
  method <- match.arg(method)
  type <- match.arg(type)
  stationary.test <- switch(method,adf = adf.test(x,nlag,output),
                            pp = pp.test(x,type, lag.short,output),
                            kpss = kpss.test(x,lag.short,output))
}