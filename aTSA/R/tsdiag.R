#' Diagnostics for ARIMA fits
#' @description Performs diagnostics for ARIMA model fitted by \code{\link{arima}} or 
#' \code{\link{estimate}} with output of diagnostic plots.
#' @param object the result of an \code{arima} or \code{estimate} fit.
#' @param lag.seq the sequence of lag to calculate the Ljung-Box test statistics. The default 
#' is \code{NULL}.
#' @details This function is similar to \code{\link{ts.diag}} in \code{stats} package, but with 
#' one more diagnostic plot for the normality of residuals. Also, the default sequence of lags  
#' for a Ljung-Box test is set to be \code{seq(4,24,by = 4)} if sample size \eqn{n > 24}, 
#' otherwise \code{seq(1,n,4)}. This function has been automatically implemented in 
#' \code{\link{estimate}} function.
#' 
#' Diagnostics are plotted, including the ACF plot, PACF plot, p.value of 
#' white noise checking plot, and Q-Q plot for residuals.
#' @return A matrix for the result of white noise checking by Ljung-Box test. 
#' @author Debin Qiu
#' @examples x <- arima.sim(list(order = c(3,0,0),ar = c(0.2,0.4,-0.15)),n = 100)
#' fit <- estimate(x,p = 3) # same as fit <- arima(x,order = c(3,0,0))
#' ts.diag(fit)
#' @importFrom stats residuals
#' @importFrom stats Box.test
#' @importFrom stats acf
#' @importFrom stats pacf
#' @importFrom stats qqnorm
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics mtext
#' @importFrom stats qqline
#' @export 
ts.diag <- function(object,lag.seq = NULL)
{
  if (class(object) != "Arima" && class(object) != "estimate")
    stop("'object' should be 'Arima' or 'estimate' class estimated
         by arima() or estimate()")
  class(object) <- "Arima"
  res <- residuals(object)
  n <- length(res)
  if (is.null(lag.seq))
    lag.seq <- if (n < 24) seq(1,n,by = 4) else seq(4,24,by = 4)
  QB <- p.value <- numeric(length(lag.seq))
  for (i in 1:length(lag.seq)) {
    k <- lag.seq[i]
    LB.test <- Box.test(res,lag = k)
    QB[i] <- LB.test$statistic
    p.value[i] <- LB.test$p.value
  }
  WNcheck <- matrix(c(lag.seq,QB,p.value),ncol = 3)
  colnames(WNcheck) <- c("lag","LB","p.value")
  op <- par(mfrow = c(2,2),oma=c(0,0,3,0))
  acf(res,main = "")
  pacf(res,main = "")
  plot(lag.seq,p.value,ylim = c(0,1),xlab = "Lag",ylab = "WN Prob.")
  abline(h = 0.05,col = 2,lty = 2)
  qqnorm(res,main = ""); qqline(res,col = 2)
  par(op)
  mtext(paste("Residual Diagnostics Plots"))
  ts.diag <- WNcheck
}