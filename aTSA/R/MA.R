#' Moving Average Filter
#' @description Applies moving average filter to estimate the linear trend or 
#' nonseasonal pattern.
#' @param x a numeric vector or univariate time series.
#' @param nlag the number of period to calculate the average. The default is \code{NULL}.
#' @param plot a logical value indicating to print out the plot. The default is \code{TRUE}.
#' @details The moving average filter uses the unweight mean of (2*\code{nlag} + 1) adjacent 
#' observations. That is, 
#' \deqn{hat{X}[t] = (X[t - nlag] + ... + X[t] + ...+ X[t + nlag])/(2*nlag + 1)}
#' for \eqn{nlag < t < n - nlag}.
#' For the values at the boundary \eqn{t \le nlag} or \eqn{n - nlag \le t \le n}, you can 
#' refer to Equation (7) in Qiu \emph{et al.,} (2013) for details of calculations.
#' The default method for choosing the optimal \code{nlag} uses the rule-of-thumb 
#' criterion proposed by Qiu, \emph{et al}., (2013), in which they showed that the moving 
#' average
#' is a special case of local linear estimator in the sense that the kernel function is the
#' uniform one, and the moving average period \code{nlag} is a function of bandwidth. Thus,
#' choosing the optimal \code{nlag} is equivalent to choosing the optimal bandwidth in local
#' linear regression.  
#' 
#' The plot of original values v.s fitted values will be displayed if \code{plot = TRUE}.
#' @return A list with class "\code{MA}" containing the following components:
#' \item{estimate}{the smoothed values.}
#' \item{nlag}{the period used to compute the average.}
#' \item{accurate}{the accurate measurements.}
#' 
#' @author Debin Qiu
#' @references
#' D. Qiu, Q. Shao, and L. Yang (2013), Efficient inference for autoregressive 
#' coefficient in the presence of trend. \emph{Journal of Multivariate Analysis}
#' 114, 40-53.
#'
#' P.J. Brockwell, R.A. Davis, Time Series: Theory and Methods, second ed., 
#' Springer, New York, 1991.
#' @examples x <- arima.sim(list(order = c(1,0,0),ar = 0.4),n = 100)
#' y <- 5*(1:100)/100 + x
#' MA(y)
#'
#'# moving average filter for co2 data 
#' MA(co2)
#' @importFrom stats lm
#' @importFrom stats is.ts
#' @importFrom stats as.ts
#' @importFrom stats ts
#' @importFrom stats resid
#' @importFrom stats time
#' @importFrom stats frequency
#' @importFrom stats var
#' @importFrom graphics plot
#' @importFrom graphics par
#' @importFrom graphics lines
#' @importFrom graphics legend
#' @export
MA <- function(x,nlag = NULL,plot = TRUE)
{
  DNAME <- deparse(substitute(x))
  if (!is.ts(x))
    x <- as.ts(x[is.finite(x)])
  n <- length(x)
  if (n < 1L)
    stop("invalid length of 'x'")
  if (is.null(nlag)) {
    xn <- (1:n)/n
    lm.fit <- lm(x ~ 1 + xn + I(xn^2) + I(xn^3))
    denominator <- 4 * (lm.fit[[1]][3])^2 + 12 * lm.fit[[1]][3] * 
      lm.fit[[1]][4] + 12 * (lm.fit[[1]][4])^2
    q.n <- floor((n)^(4/5) * (9/2)^(1/5) * ((var(resid(lm.fit))/denominator)^(1/5)))
    nlag <- ifelse(q.n < n, min(q.n, n - q.n), floor(n^(4/5)/2))
  }
  else {
    if (nlag < 1L || nlag > n || nlag%%1 != 0)
      stop("'nlag' must be an integer between [1,n]")
  }
  x.hat <- numeric(n)
  for (i in 1:n) {
    if (i <= nlag) x.hat[i] <- mean(x[1:(i + nlag)]) 
    else if (i > nlag && i < (n - nlag)) x.hat[i] <- mean(x[(i - nlag):(i + nlag)])
    else x.hat[i] <- mean(x[(i - nlag):n])
  }
  names(nlag) <- NULL
  x.hat <- ts(x.hat,start = time(x)[1],frequency = frequency(x))
  result <- list(estimate = x.hat, nlag = nlag, 
                 accurate = accurate(x,x.hat,1,output = FALSE)) 
  if (plot) {
    par(mar = c(5, 4, 1.4, 0.2))
    plot(x,xlab = "time",ylab = DNAME)
    lines(x.hat,col = 2)
    legend('topright',legend = c("original","fitted"),
           col = c(1,2),lty = c(1,1),bty ="n",horiz=TRUE)
  }
  class(result) <- "MA"
  return(result)
}