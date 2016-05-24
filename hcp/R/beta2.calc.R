beta2.calc <-
function(x, y, n, j, k, e1, e2)
{
	aa <- wmat2(x, y, n, j, k, e1, e2)
	W <- aa$w
	bb <- rvec2(x, y, n, j, k, e1, e2)
	R <- bb$r
	beta <- solve(W, R)
	list(B = beta)
}
