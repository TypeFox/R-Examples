d0 <-
function(z, y) {
	if(length(y) != length(z)) y <- rep(y, length(z))
	pr <- y > 0 & y < 1
	npr <- !pr
	ans <- vector('numeric', length(y))
	y[pr] <- qnorm(y[pr])
	y[npr] <- 2 * y[npr] - 1
	ans[pr] <- dnorm(y[pr], z[pr], log = TRUE)
	ans[npr] <- pnorm(y[npr] * z[npr], log.p = TRUE)
	ans
}
