Box.Ljung.Test = function(z, lag = NULL, main = NULL){
  if(is.null(lag)){lag = 10}
  k = lag
  n = length(z)
  aux = acf(z, plot = FALSE, lag.max = k, na.action = na.pass)
  p.value = vector("numeric")
  Q = vector("numeric")
  for(j in 1:k){
    rho = aux$acf[2:(j+1),,1]
    Q[j] = sum(n*(n+2)*rho^2/(n-1:j))
    p.value[j] = 1-pchisq(Q[j], df = j)
  }
  if(is.null(main)){main = expression("p values for Ljung-Box statistic")}
  plot(p.value ~ c(1:k), ylim = c(0,1), bty = "n", las = 1, lwd = 2, xlim = c(0,k), main = main, xlab = "Lag", ylab = "p-value", pch = 20)
  abline(h = 0.05, lty = 2, col = "blue")
}
