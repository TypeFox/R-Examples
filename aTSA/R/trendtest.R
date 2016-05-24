#' Trend Test
#' @description Performs an approximate Cox-Stuart or Difference-Sign trend test.
#' @param x a numeric vector or univariate time series.
#' @param method test method. The default is \code{method = "cox.stuart"}.
#' @param plot a logical value indicating to display the plot of data. 
#' The default is \code{FALSE}.
#' @details Cox-Stuart or Difference-Sign test is used to test whether the data have a
#' increasing or decreasing trend. They are useful to detect the linear or nonlinear trend.
#' The Cox-Stuart test is constructed as follows. 
#' For the given data \eqn{x[1],...,x[t]}, one can divide them into two sequences with 
#' equal number of observations cutted in the midpoint and then take the paired difference, 
#' i.e., \eqn{D = x[i] - x[i+c], i = 1, ..., floor(n/2)}, where \eqn{c} is the index of 
#' midpoint. Let \eqn{S} be the number of positive or negative values in \eqn{D}. Under the
#' null hypothesis that data have no trend, for large \eqn{n} = length(x), \eqn{S} is 
#' approximately distributed as \eqn{N(n/2,n/4)}, such that one can immediately obtain 
#' the p value. The exact Cox-Stuart trend test can be seen in \code{\link[snpar]{cs.test}} of
#' \code{snpar} package.
#' 
#' The Difference-Sign test is constructed as the similar way as Cox-Stuart test. We first 
#' let \eqn{D = x[i] - x[i - 1]} for \eqn{i = 2, ..., n} and then
#' count the number of positive or negative values in \eqn{D}, defined as \eqn{S}. 
#' Under the null hypothesis, \eqn{S} is approximately distributed as 
#' \eqn{N((n-1)/2,(n+1)/12)}. Thus, p-value can be calculated based on the null distribution.
#' 
#' @note Missing values are removed.
#' @return A list with class "\code{htest}" containing the following components:
#' \item{data.name}{a character string giving the names of the data.}
#' \item{method}{the type of test applied.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{p.value}{the p-value for the test.}
#' \item{statistic}{the value of the test statistic with a name describing it.}
#' @author Debin Qiu
#' @references 
#' D.R. Cox and A. Stuart (1955). Some quick sign tests for trend in location 
#' and dispersion. \emph{Biometrika}, Vol. 42, pp. 80-95.
#' 
#' P.J. Brockwell, R.A. Davis, Time Series: Theory and Methods, second ed., 
#' Springer, New York, 1991. (p. 37)
#' 
#' @examples x <- rnorm(100)
#' trend.test(x,plot = TRUE) # no trend
#' 
#' x <- 5*(1:100)/100
#' x <- x + arima.sim(list(order = c(1,0,0),ar = 0.4),n = 100)
#' trend.test(x,plot = TRUE) # increasing trend
#' @importFrom stats pnorm
#' @importFrom stats plot.ts
#' @importFrom stats setNames
#' @export 
trend.test <- function(x,method = c("cox.stuart","diff.sign"),plot = FALSE)
{
  DNAME <- deparse(substitute(x))
  if (NCOL(x) > 1L)
    stop("'x' must be a numeric vector or univariate time series")
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L)
    stop("not enough 'x' data")
  if (plot)
    plot.ts(x,main = paste("Plot of",DNAME))
  method <- match.arg(method)
  METHOD <- switch(method,cox.stuart = "Approximate Cox-Stuart trend test",
                   diff.sign = "Approximate Difference-Sign trend test")
  c <- as.integer(ifelse(n%%2 == 0, n/2, (n - 1)/2))
  Diff <- switch(method,diff.sign = x[-1] - x[-n],
              cox.stuart = if (n%%2 == 0)  x[1:c] - x[(c + 1):n] 
              else x[1:c] - x[(c + 2):n])
  n <- length(Diff[Diff != 0])
  D.pos <- sum(Diff > 0)
  D.neg <- sum(Diff < 0)
  D <- max(D.pos, D.neg)
  Z <- switch(method,diff.sign = (D - (n - 1)/2)/sqrt((n + 1)/12),
              cox.stuart = (D - n/2)/sqrt(n/4))
  p.val <- min(pnorm(Z),1 - pnorm(Z))
  ALTERNATIVE <- if (D.pos > D.neg) "data have a increasing trend" 
                 else "data have a decreasing trend"
  STATISTIC <- setNames(D,ifelse(D.pos > D.neg,"D+","D-"))
  RVAL <- list(data.name = DNAME, method = METHOD, alternative = ALTERNATIVE, 
               p.value = p.val, statistic = STATISTIC)
  class(RVAL) <- "htest"
  return(RVAL)
}