plot.ml_g_fit <- function(x, ...) {
  e.hat <- residuals(x)
  e.hat.ss <- residuals(x, type = "ss")
  y.hat <- fitted(x)
  n <- nrow(x$X)
  par(mfrow = c(2,2), mar = c(4,4,2,1), las = 1)
  scatter.smooth(x = y.hat, y = x$y,         ### Plot 1
                 ylab = "Observations", xlab = "Fitted Values")
  abline(a = 0, b = 1)
  scatter.smooth(x = y.hat, y = e.hat,       ### Plot 2
                 ylab = "Residuals", xlab = "Fitted Values") 
  abline(h = 0)
  scatter.smooth(x = y.hat, y = abs(sqrt(e.hat.ss)), ### Plot 3
                 ylab = "Sqrt (Abs( Stand. Res.))",
                 xlab = "Fitted Values") 
  plot(x = sort(e.hat.ss), y = qnorm((1:n)/(n+1)), ### Plot 4
       xlab = "Standardized Studentized Residuals",
       ylab = "Normal Quantiles")
}
