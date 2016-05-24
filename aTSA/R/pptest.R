#' Phillips-Perron Test
#' @description Performs the Phillips-Perron test for the null hypothesis of a unit root of 
#' a univariate time series \code{x} (equivalently, \code{x} is a non-stationary time series).
#' @param x a numeric vector or univariate time series.
#' @param type the type of Phillips-Perron test. The default is \code{Z_rho}.
#' @param lag.short a logical value indicating whether the parameter of lag to calculate the 
#' statistic is a short or long term. The default is a short term.
#' @param output a logical value indicating to print the results in R console. The default
#' is \code{TRUE}.
#' @details Compared with the Augmented Dickey-Fuller test, Phillips-Perron test makes 
#' correction to the test statistics and is robust to the unspecified autocorrelation 
#' and heteroscedasticity in the errors. There are two types of test statistics, 
#' \eqn{Z_{\rho}} and \eqn{Z_{\tau}}, which have the same asymptotic distributions
#' as Augmented Dickey-Fuller test statistic, \code{ADF}. The calculations of each type 
#' of the Phillips-Perron test can be see in the reference below. If the 
#' \code{lag.short = TRUE}, we use the default number of Newey-West lags 
#' \eqn{floor(4*(length(x)/100)^0.25)}, 
#' otherwise \eqn{floor(12*(length(x)/100)^0.25)} to calculate the test statistics. 
#' In order to calculate the test statistic, we consider
#' three types of linear regression models. The first type (\code{type1}) is the one
#'  with no drift and linear trend with respect to time:
#' \deqn{x[t] = \rho*x[t-1] + e[t],}
#' where \eqn{e[t]} is an error term. 
#' The second type (\code{type2}) is the one with drift but no linear trend:
#' \deqn{x[t] = \mu + \rho*x[t-1] + e[t].}
#' The third type (type3) is the one with both drift and linear trend:
#' \deqn{x[t] = \mu + \alpha*t + \rho*x[t-1] + e[t].}
#' The p.value is calculated by the interpolation of test statistics from the critical values
#' tables (Table 10.A.1 for \code{Z_rho} and 10.A.2 for \code{Z_tau} in Fuller (1996))
#' with a given sample size \eqn{n} = length(\code{x}). 
#' 
#' @note Missing values are removed.
#' @return A matrix for test results with three columns (\code{lag},\code{Z_rho} 
#' or \code{Z_tau}, \code{p.value}) and three rows (\code{type1}, \code{type2}, \code{type3}).
#' Each row is the test results (including lag parameter, test statistic and p.value) for 
#' each type of linear equation. 
#' @author Debin Qiu
#' @references 
#' Phillips, P. C. B.; Perron, P. (1988). Testing for a Unit Root in Time Series Regression. 
#' \emph{Biometrika}, 75 (2): 335-346.
#' 
#' Fuller, W. A. (1996). Introduction to statistical time series, second ed., Wiley, New York. 
#' @seealso \code{\link{adf.test}}, \code{\link{kpss.test}}, \code{\link{stationary.test}} 
#' @examples 
#' # PP test for ar(1) process
#' x <- arima.sim(list(order = c(1,0,0),ar = 0.2),n = 100)
#' pp.test(x)
#' 
#' # PP test for co2 data
#' pp.test(co2)
#' @importFrom stats embed
#' @importFrom stats residuals
#' @importFrom stats approx
#' @importFrom stats lm
#' @export
pp.test <- function(x,type = c("Z_rho","Z_tau"),
                    lag.short = TRUE,output = TRUE)
{
  type <- match.arg(type)
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
  lag <- ifelse(lag.short,4,12)
  q <- floor(lag*(n/100)^0.25)  
  stat <- function(model,m) {
    index <- ifelse (m > 1, 2, 1)
    res1 <- residuals(model)
    rho.hat <- summary(model)$coefficients[index,1]
    sig.hat <- summary(model)$coefficients[index,2]
    s2 <- sum(res1^2)/(n - m)
    gamma <- numeric(q + 1)
    for (i in 1:(q + 1)) {
      u <- embed(res1,i)
      gamma[i] = sum(u[,1]*u[,i])/n
    }
    lambda2 <- gamma[1] + 2*sum((1 - 1:q/(q + 1))*gamma[-1])
    Z_r <- n*(rho.hat - 1) - n^2*sig.hat^2/s2*(lambda2 - gamma[1])/2
    Z_t <- sqrt(gamma[1]/lambda2)*(rho.hat - 1)/sig.hat - 
           (lambda2 - gamma[1])*n*sig.hat/(2*sqrt(lambda2*s2))
    return(c(Z_r,Z_t))
  }
  pvalue <- function(table,stat) {
    Ncol <- ncol(table)
    Size <- c(25, 50, 100, 250, 500, 1000)
    Percnt <- c(.01, .025, .05, .10, .50, .90, .95, .975, .99)
    intplSize <- numeric(Ncol)
    for (j in 1:Ncol) intplSize[j] <- approx(Size, table[,j], n, rule = 2)$y
    approx(intplSize, Percnt, stat, rule = 2)$y
  } 
  table_r1 <- rbind(c(-11.8,-9.3,-7.3,-5.3,-0.82,1.01,1.41,1.78,2.28),
                    c(-12.8,-9.9,-7.7,-5.5,-0.84,0.97,1.34,1.69,2.16),
                    c(-13.3,-10.2,-7.9,-5.6,-0.85,0.95,1.31,1.65,2.09),
                    c(-13.6,-10.4,-8.0,-5.7,-0.86,0.94,1.29,1.62,2.05),
                    c(-13.7,-10.4,-8.0,-5.7,-0.86,0.93,1.29,1.61,2.04),
                    c(-13.7,-10.5,-8.1,-5.7,-0.86,0.93,1.28,1.60,2.03))
  table_r2 <- rbind(c(-17.2,-14.6,-12.5,-10.2,-4.22,-0.76,0.00,0.64,1.39),
                    c(-18.9,-15.7,-13.3,-10.7,-4.29,-0.81,-0.07,0.53,1.22),
                    c(-19.8,-16.3,-13.7,-11.0,-4.32,-0.83,-0.11,0.47,1.13),
                    c(-20.3,-16.7,-13.9,-11.1,-4.34,-0.84,-0.13,0.44,1.08),
                    c(-20.5,-16.8,-14.0,-11.2,-4.35,-0.85,-0.14,0.42,1.07),
                    c(-20.6,-16.9,-14.1,-11.3,-4.36,-0.85,-0.14,0.41,1.05))
  table_r3 <- rbind(c(-22.5,-20.0,-17.9,-15.6,-8.49,-3.65,-2.51,-1.53,-0.46),
                    c(-25.8,-22.4,-19.7,-16.8,-8.80,-3.71,-2.60,-1.67,-0.67),
                    c(-27.4,-23.7,-20.6,-17.5,-8.96,-3.74,-2.63,-1.74,-0.76),
                    c(-28.5,-24.4,-21.3,-17.9,-9.05,-3.76,-2.65,-1.79,-0.83),
                    c(-28.9,-24.7,-21.5,-18.1,-9.08,-3.76,-2.66,-1.80,-0.86),
                    c(-29.4,-25.0,-21.7,-18.3,-9.11,-3.77,-2.67,-1.81,-0.88))
  table_t1 <- rbind(c(-2.65,-2.26,-1.95,-1.60,-0.47,0.92,1.33,1.70,2.15),
                  c(-2.62,-2.25,-1.95,-1.61,-0.49,0.91,1.31,1.66,2.08),
                  c(-2.60,-2.24,-1.95,-1.61,-0.50,0.90,1.29,1.64,2.04),
                  c(-2.58,-2.24,-1.95,-1.62,-0.50,0.89,1.28,1.63,2.02),
                  c(-2.58,-2.23,-1.95,-1.62,-0.50,0.89,1.28,1.62,2.01),
                  c(-2.58,-2.23,-1.95,-1.62,-0.51,0.89,1.28,1.62,2.01))
  table_t2 <- rbind(c(-3.75,-3.33,-2.99,-2.64,-1.53,-0.37,0.00,0.34,0.71),
                  c(-3.59,-3.23,-2.93,-2.60,-1.55,-0.41,-0.04,0.28,0.66),
                  c(-3.50,-3.17,-2.90,-2.59,-1.56,-0.42,-0.06,0.26,0.63),
                  c(-3.45,-3.14,-2.88,-2.58,-1.56,-0.42,-0.07,0.24,0.62),
                  c(-3.44,-3.13,-2.87,-2.57,-1.57,-0.44,-0.07,0.24,0.61),
                  c(-3.42,-3.12,-2.86,-2.57,-1.57,-0.44,-0.08,0.23,0.60))
  table_t3 <- rbind(c(-4.38,-3.95,-3.60,-3.24,-2.14,-1.14,-0.81,-0.50,-0.15),
                  c(-4.16,-3.80,-3.50,-3.18,-2.16,-1.19,-0.87,-0.58,-0.24),
                  c(-4.05,-3.73,-3.45,-3.15,-2.17,-1.22,-0.90,-0.62,-0.28),
                  c(-3.98,-3.69,-3.42,-3.13,-2.18,-1.23,-0.92,-0.64,-0.31),
                  c(-3.97,-3.67,-3.42,-3.13,-2.18,-1.24,-0.93,-0.65,-0.32),
                  c(-3.96,-3.67,-3.41,-3.13,-2.18,-1.25,-0.94,-0.66,-0.32))
 m1 <- lm(yt ~ yt1 - 1)
 m2 <- lm(yt ~ yt1)
 m3 <- lm(yt ~ yt1 + t)
 STAT <- rbind(stat(m1,1),stat(m2,2),stat(m3,3))
 PVAL.r <- c(pvalue(table_r1,STAT[1,1]),pvalue(table_r2,STAT[2,1]),
              pvalue(table_r3,STAT[3,1]))
 PVAL.t <- c(pvalue(table_t1,STAT[1,2]),pvalue(table_t2,STAT[2,2]),
             pvalue(table_t3,STAT[3,2])) 
 STATISTIC <- switch(type, Z_rho = STAT[,1], Z_tau = STAT[,2])
 PVALUE <- switch(type, Z_rho = PVAL.r, Z_tau = PVAL.t)
 if (output) {
   NAME <- switch(type, Z_rho = list(c(""),c("lag","Z_rho","p.value")),
                  Z_tau = list(c(""),c("lag","Z_tau","p.value")))
   cat("Phillips-Perron Unit Root Test","\n")
   cat("alternative: stationary","\n","\n")
   cat("Type 1: no drift no trend","\n")
   print(matrix(c(q,STATISTIC[1],PVALUE[1]),1,3,dimnames = NAME), digits = 3)
   cat("-----","\n","Type 2: with drift no trend","\n")
   print(matrix(c(q,STATISTIC[2],PVALUE[2]),1,3,dimnames = NAME), digits = 3)
   cat("-----","\n","Type 3: with drift and trend","\n")
   print(matrix(c(q,STATISTIC[3],PVALUE[3]),1,3,dimnames = NAME), digits = 3)
   cat("---------------","\n")
   cat("Note: p-value = 0.01 means p.value <= 0.01","\n")
 }
 coln <- switch(type,Z_rho = c("lag","Z_rho","p.value"),
                Z_tau = c("lag","Z_tau","p.value"))
 pp.test <- matrix(c(rep(q,3),STATISTIC,PVALUE),3,3,
                  dimnames = list(paste("type",1:3),coln))                                 
}