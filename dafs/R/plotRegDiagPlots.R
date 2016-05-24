plotRegDiagPlots = function(fit){
  par(mfrow = c(2,2), mar = c(3,3,2,1), mgp = c(1,0,0), oma = c(0,0,0,0))
  res = residuals(fit)
  pred = fitted(fit)
  plot(pred, res, axes = FALSE, xlab = "Fitted values", ylab  = "Residuals")
  abline(h = 0, col = "grey", lty = 2)
  box()

  mx = mean(res)
  sx = sd(res)
  hist(res, prob = TRUE, axes = FALSE, main = "", xlab = "Residuals",
       ylab = "Density")
  x = seq(min(res), max(res), length = 200)
  y = dnorm(x, mx, sx)
  lines(x, y, lty = 2, col = "grey")
  box()

  n = length(res)
  z = qnorm(ppoints(n))
  res = sort(res)
  fit = lm(res~z)
  qqnorm(res, axes = FALSE, main = "")
  abline(fit)
  fit = lm(res~z)
  abline(fit)
  box()
}
