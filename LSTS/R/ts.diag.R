ts.diag = function(x, lag = 10, cex = 0.5){
  Z = (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
  op = par(mfrow = c(3,1), cex = cex)
  plot(Z ~ time(Z), ylim = c(-3.5,3.5), type = "h", lwd = 1, main = expression("Standardized Residuals"), xlab = "", ylab = "", bty = "n", las = 1)
  abline(h = 0)
  abline(h = +2, col = 4, lwd = 1, lty = 2)
  abline(h = -2, col = 4, lwd = 1, lty = 2)
  abline(h = +3, col = 4, lwd = 1, lty = 3)
  abline(h = -3, col = 4, lwd = 1, lty = 3)
  acf(x, lag.max = lag, main = expression("ACF of Residuals"), lwd = 1, las = 1, col = 1, bty = "n", na.action = na.pass)
  Box.Ljung.Test(x, lag = lag)
  par(op)
}