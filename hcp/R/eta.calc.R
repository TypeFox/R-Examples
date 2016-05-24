eta.calc <-
function(x, y, n, j, k, theta)
{
	jp1 <- j + 1
	kp1 <- k + 1
	a0 <- theta[1]
	a1 <- theta[2]
	b1 <- theta[3]
	c1 <- theta[4]
	b0 <- a0 + (a1 - b1) * x[j]
	c0 <- b0 + (b1 - c1) * x[k]
	rss1 <- sum((y[1:j] - a0 - a1 * x[1:j])^2)
	rss2 <- sum((y[jp1:k] - b0 - b1 * x[jp1:k])^2)
	rss3 <- sum((y[kp1:n] - c0 - c1 * x[kp1:n])^2)
	e1 <- (n - k + j)/(rss1 + rss3)
	e2 <- (k - j)/rss2
	list(eta1 = e1, eta2 = e2)
}
