d2 <-
function(z, y) {
	pr <- y > 0 & y < 1
	npr <- !pr
	ans <- vector('numeric', length(y))
	y[npr] <- 2 * y[npr] - 1
	ans[pr] <- -1
	a <- dnorm(z[npr]) ^ 2 / pnorm(y[npr] * z[npr]) ^ 2
	b <- y[npr] * z[npr] * dnorm(z[npr]) / pnorm(y[npr] * z[npr])
	ans[npr] <- -a - b
	ans
}
