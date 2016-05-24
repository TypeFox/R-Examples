d1 <-
function(z, y) {
	pr <- y > 0 & y < 1
	npr <- !pr
	ans <- vector('numeric', length(y))
	y[pr] <- qnorm(y[pr])
	y[npr] <- 2 * y[npr] - 1
	ans[pr] <- y[pr] - z[pr]
	ans[npr] <- y[npr] * dnorm(z[npr]) / pnorm(y[npr] * z[npr])
	ans
}

