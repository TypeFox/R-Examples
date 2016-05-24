#' Augmented Dickey-Fuller Test
#' @description Performs the Augmented Dickey-Fuller test for the null hypothesis 
#' of a unit root of a univarate time series \code{x} (equivalently, \code{x} is a 
#' non-stationary time series).
#' @param x a numeric vector or univariate time series.
#' @param nlag the lag order with default to calculate the test statistic. See details for 
#' the default.
#' @param output a logical value indicating to print the test results in R console. 
#' The default is \code{TRUE}.
#' 
#' @details The Augmented Dickey-Fuller test incorporates 
#' three types of linear regression models. The first type (\code{type1}) is a linear model 
#' with no drift and linear trend with respect to time:
#' \deqn{dx[t] = \rho*x[t-1] + \beta[1]*dx[t-1] + ... + \beta[nlag - 1]*dx[t - nlag + 1]
#' +e[t],}
#' where \eqn{d} is an operator of first order difference, i.e., 
#' \eqn{dx[t] = x[t] - x[t-1]}, and \eqn{e[t]} is an error term.
#' 
#' The second type (\code{type2}) is a linear model with drift but no linear trend:
#' \deqn{dx[t] = \mu  + \rho*x[t-1] + \beta[1]*dx[t-1] + ... + 
#' \beta[nlag - 1]*dx[t - nlag + 1] +e[t].}
#' 
#' The third type (\code{type3}) is a linear model with both drift and linear trend:
#' \deqn{dx[t] = \mu + \beta*t + \rho*x[t-1] + \beta[1]*dx[t-1] + ... + 
#' \beta[nlag - 1]*dx[t - nlag + 1] +e[t].}
#' 
#' We use the default \code{nlag = floor(4*(length(x)/100)^(2/9))} to 
#' calcuate the test statistic. 
#' The Augmented Dickey-Fuller test statistic is defined as 
#' \deqn{ADF = \rho.hat/S.E(\rho.hat),}
#' where \eqn{\rho.hat} is the coefficient estimation
#' and \eqn{S.E(\rho.hat)} is its corresponding estimation of standard error for each 
#' type of linear model. The p.value is 
#' calculated by interpolating the test statistics from the corresponding critical values 
#' tables (see Table 10.A.2 in Fuller (1996)) for each type of linear models with given 
#' sample size \eqn{n} = length(\code{x}). 
#' The Dickey-Fuller test is a special case of Augmented Dickey-Fuller test 
#' when \code{nlag} = 2. 
#'  
#' @note Missing values are removed.
#' 
#' @return A list containing the following components: 
#'    \item{type1}{a matrix with three columns: \code{lag}, \code{ADF}, \code{p.value}, 
#'    where \code{ADF} is the Augmented Dickey-Fuller test statistic.}
#'    \item{type2}{same as above for the second type of linear model.}
#'    \item{type3}{same as above for the third type of linear model.}
#' 
#'  @author Debin Qiu 
#'  @references 
#'  Fuller, W. A. (1996). Introduction to Statistical Time Series, second ed., New York: 
#'  John Wiley and Sons.
#'  @seealso \code{\link{pp.test}}, \code{\link{kpss.test}}, \code{\link{stationary.test}}
#'  @examples
#'  # ADF test for AR(1) process
#'x <- arima.sim(list(order = c(1,0,0),ar = 0.2),n = 100)
#'adf.test(x)
#'# ADF test for co2 data
#'adf.test(co2)
#'@importFrom stats embed
#'@importFrom stats lm
#'@importFrom stats update
#'@importFrom stats approx
#'@export
adf.test <- function(x, nlag = NULL, output = TRUE)
{
  if (NCOL(x) > 1)
    stop("'x' must be a numeric vector or univariate time series")
  x <- x[is.finite(x)]
  if (length(x) < 1L)
    stop("invalid length of 'x'") 
  if (!is.null(nlag) && (nlag%%1 != 0 || nlag < 0))
    stop("'nlag' must be a positive integer number")
  nlag <- ifelse(is.null(nlag),floor(4*(length(x)/100)^(2/9)) + 1,nlag)
  y <- diff(x)
  n <- length(y)
  result1 <- result2 <- result3 <- NULL
  for (i in 1:nlag) {
    z <- embed(y,i)
    yt <- z[,1]
    t <- i:n
    xt1 <- x[t]
    m1 <- lm(yt ~ xt1 - 1)
    m2 <- lm(yt ~ xt1)
    m3 <- lm(yt ~ xt1 + t)
    if (i > 1) {
      yt1 <- z[,2:i]
      m1 <- update(m1,.~. + yt1)
      m2 <- update(m2,.~. + yt1)
      m3 <- update(m3,.~. + yt1)
    }
    STAT <- c(summary(m1)$coefficients[1,1]/summary(m1)$coefficients[1,2],
              summary(m2)$coefficients[2,1]/summary(m2)$coefficients[2,2],
              summary(m3)$coefficients[2,1]/summary(m3)$coefficients[2,2]) 
    table1 <- rbind(c(-2.65,-2.26,-1.95,-1.60,-0.47,0.92,1.33,1.70,2.15),
                    c(-2.62,-2.25,-1.95,-1.61,-0.49,0.91,1.31,1.66,2.08),
                    c(-2.60,-2.24,-1.95,-1.61,-0.50,0.90,1.29,1.64,2.04),
                    c(-2.58,-2.24,-1.95,-1.62,-0.50,0.89,1.28,1.63,2.02),
                    c(-2.58,-2.23,-1.95,-1.62,-0.50,0.89,1.28,1.62,2.01),
                    c(-2.58,-2.23,-1.95,-1.62,-0.51,0.89,1.28,1.62,2.01))
    table2 <- rbind(c(-3.75,-3.33,-2.99,-2.64,-1.53,-0.37,0.00,0.34,0.71),
                    c(-3.59,-3.23,-2.93,-2.60,-1.55,-0.41,-0.04,0.28,0.66),
                    c(-3.50,-3.17,-2.90,-2.59,-1.56,-0.42,-0.06,0.26,0.63),
                    c(-3.45,-3.14,-2.88,-2.58,-1.56,-0.42,-0.07,0.24,0.62),
                    c(-3.44,-3.13,-2.87,-2.57,-1.57,-0.44,-0.07,0.24,0.61),
                    c(-3.42,-3.12,-2.86,-2.57,-1.57,-0.44,-0.08,0.23,0.60))
    table3 <- rbind(c(-4.38,-3.95,-3.60,-3.24,-2.14,-1.14,-0.81,-0.50,-0.15),
                    c(-4.16,-3.80,-3.50,-3.18,-2.16,-1.19,-0.87,-0.58,-0.24),
                    c(-4.05,-3.73,-3.45,-3.15,-2.17,-1.22,-0.90,-0.62,-0.28),
                    c(-3.98,-3.69,-3.42,-3.13,-2.18,-1.23,-0.92,-0.64,-0.31),
                    c(-3.97,-3.67,-3.42,-3.13,-2.18,-1.24,-0.93,-0.65,-0.32),
                    c(-3.96,-3.67,-3.41,-3.13,-2.18,-1.25,-0.94,-0.66,-0.32))
    pvalue <- function(table,stat) {
      Ncol <- ncol(table)
      Size <- c(25, 50, 100, 250, 500, 100000)
      Percnt <- c(.01, .025, .05, .10, .50, .90, .95, .975, .99)
      intplSize <- numeric(Ncol)
      for (j in 1:Ncol) intplSize[j] <- approx(Size, table[,j], n, rule = 2)$y
      approx(intplSize, Percnt, stat, rule = 2)$y
    }
    PVAL <- c(pvalue(table1,STAT[1]),pvalue(table2,STAT[2]),pvalue(table3,STAT[3]))
    result1 <- rbind(result1,c(i - 1,STAT[1],PVAL[1]))
    result2 <- rbind(result2,c(i - 1,STAT[2],PVAL[2]))
    result3 <- rbind(result3,c(i - 1,STAT[3],PVAL[3]))
  }
 colnames(result1) <- colnames(result2) <- colnames(result3) <- c("lag","ADF","p.value")
 if (output) {
   cat("Augmented Dickey-Fuller Test","\n")
   cat("alternative: stationary","\n","\n")
   cat("Type 1: no drift no trend","\n")
   print(result1,digits = 3)
   cat("Type 2: with drift no trend","\n")
   print(result2,digits = 3)
   cat("Type 3: with drift and trend","\n")
   print(result3,digits = 3)
   cat("----","\n")
   cat("Note: in fact, p.value = 0.01 means p.value <= 0.01","\n")
 } 
 adf.test <- list(type1 = result1,type2 = result2,type3 = result3)
}