#' Kwiatkowski-Phillips-Schmidt-Shin Test
#' @description Performs Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test for the null 
#' hypothesis that \code{x} is a stationary univariate time series.
#' @param x a numeric vector or univariate time series.
#' @param lag.short a logical value indicating whether the parameter of lag to calculate 
#' the test statistic is a short or long term. The default is a short term. See details.
#' @param output a logical value indicating to print out the results in R console. 
#' The default is \code{TRUE}.
#' @details The Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test tends to decompose the time 
#' series into the sum of a deterministic trend, a random walk, and a stationary error: 
#' \deqn{x[t] = \alpha*t + u[t] + e[t],}
#' where \eqn{u[t]} satisfies \eqn{u[t] = u[t-1] + a[t]}, and \eqn{a[t]} are i.i.d 
#' \eqn{(0,\sigma^2)}. The null hypothesis is that \eqn{\sigma^2 = 0}, which implies 
#' \code{x} is a stationary time series. In order to calculate the test statistic, 
#' we consider three types of linear regression models. 
#' The first type (\code{type1}) is the one with no drift and deterministic trend, 
#' defined as \deqn{x[t] = u[t] + e[t].}
#' The second type (\code{type2}) is the one with drift but no trend: 
#' \deqn{x[t] = \mu + u[t] + e[t].}
#' The third type (\code{type3}) is the one with both drift and trend:
#' \deqn{x[t] = \mu + \alpha*t + u[t] + e[t].}
#' The details of calculation of test statistic (\code{kpss}) can be seen in the references 
#' below. The default parameter of lag to calculate the test statistic is 
#' \eqn{max(1,floor(3*sqrt(n)/13)} for short term effect, otherwise, 
#' \eqn{max(1,floor(10*sqrt(n)/13)} for long term effect. 
#' The p.value is calculated by the interpolation of test statistic from tables of 
#' critical values (Table 5, Hobijn B., Franses PH. and Ooms M (2004)) for a given 
#' sample size \eqn{n} = length(\code{x}).
#' 
#' @note Missing values are removed.
#' @return A matrix for test results with three columns (\code{lag}, \code{kpss},  
#' \code{p.value}) and three rows (\code{type1}, \code{type2}, \code{type3}). 
#' Each row is the test results (including lag parameter, test statistic and p.value) for 
#' each type of linear regression models.
#' @author Debin Qiu
#' @references
#' Hobijn B, Franses PH and Ooms M (2004). Generalization of the KPSS-test for stationarity.
#' \emph{Statistica Neerlandica}, vol. 58, p. 482-502.
#' 
#' Kwiatkowski, D.; Phillips, P. C. B.; Schmidt, P.; Shin, Y. (1992). 
#' Testing the null hypothesis of stationarity against the alternative of a unit root.
#' \emph{Journal of Econometrics}, 54 (1-3): 159-178.
#' @seealso \code{\link{adf.test}}, \code{\link{pp.test}}, \code{\link{stationary.test}}
#' @examples 
#' # KPSS test for AR(1) process
#' x <- arima.sim(list(order = c(1,0,0),ar = 0.2),n = 100)
#' kpss.test(x)
#' 
#' # KPSS test for co2 data
#' kpss.test(co2)
#' @importFrom stats embed
#' @importFrom stats residuals
#' @importFrom stats lm
#' @importFrom stats approx
#' @export
kpss.test <- function(x, lag.short = TRUE,output = TRUE)
{
  if (NCOL(x) > 1)
    stop("'x' must be a numeric vector or univariate time series")
  x <- x[is.finite(x)]
  z <- embed(x,2)
  yt <- z[,1]
  yt1 <- z[,2]
  n <- length(yt)
  if (n < 1L)
    stop("invalid length of 'x'")
  t <- 1:n
  lag <- ifelse(lag.short,3,10)
  q <- max(1,floor(lag * sqrt(n)/13)) 
  stat <- function(model,m) {
    index <- ifelse (m > 1, 2, 1)
    res <- residuals(model)
    S2 <- cumsum(res)^2   
    gamma <- numeric(q + 1)
    for (i in 1:(q + 1)) {
      u <- embed(res,i)
      gamma[i] = sum(u[,1]*u[,i])/n
    }
    sigma2 <- gamma[1] + 2*sum((1 - 1:q/(q + 1))*gamma[-1])
    sum(S2)/(n^2*sigma2)    
  } 
  m1 <- lm(yt ~ yt1 - 1)
  m2 <- lm(yt ~ yt1)
  m3 <- lm(yt ~ yt1 + t)
  STAT <- c(stat(m1,1),stat(m2,2),stat(m3,3))
  table1 <- c(1.195,1.656,2.114,2.759)
  table2 <- c(0.347,0.463,0.574,0.739)
  table3 <- c(0.119,0.146,0.176,0.216)
  percnt <- c(.10,.05,.025,.01)
  PVAL <- c(approx(table1,percnt,STAT[1],rule = 2)$y,
            approx(table2,percnt,STAT[2],rule = 2)$y,
            approx(table3,percnt,STAT[3],rule = 2)$y)  
  if (output) {
    NAME <- list(c(""),c("lag","stat","p.value"))
    cat("KPSS Unit Root Test","\n")
    cat("alternative: nonstationary","\n","\n")
    cat("Type 1: no drift no trend","\n")
    print(matrix(c(q,STAT[1],PVAL[1]),1,3,dimnames = NAME), digits = 3) 
    cat("-----","\n","Type 2: with drift no trend","\n")
    print(matrix(c(q,STAT[2],PVAL[2]),1,3,dimnames = NAME), digits = 3)
    cat("-----","\n","Type 1: with drift and trend","\n")
    print(matrix(c(q,STAT[3],PVAL[3]),1,3,dimnames = NAME), digits = 3)
    cat("-----------","\n")
    cat("Note: p.value = 0.01 means p.value <= 0.01","\n")
    cat("    : p.value = 0.10 means p.value >= 0.10","\n")
  }
  kpss.test <- matrix(c(rep(q,3),STAT,PVAL),3,3,dimnames = 
                     list(paste("type",1:3),c("lag","kpss","p.value")))
}