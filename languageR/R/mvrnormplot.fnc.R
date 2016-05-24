`mvrnormplot.fnc` <- 
function(r = 0.9, n = 100, limits=NA) {
	require("MASS", quietly = TRUE)
  x = mvrnorm(n, mu = c(0, 0), Sigma=cbind(c(1, r), c(r, 1)))
	if (is.na(limits)) 
    plot(x, xlab = "x", ylab = "y")
	else
    plot(x, xlab = "x", ylab = "y", ylim=limits, xlim=limits)
  abline(lm(x[, 2] ~ x[, 1]))
  mtext(paste("r =", round(cor(x)[1, 2], 3)), line = 1)
	if (r > 0)
    abline(0, 1, lty = 2)
  else
    abline(0, -1, lty = 2)
}
