con.est <-
function(xj, xk, est)
{
	eta1 <- est[1]
	eta2 <- est[2]
	a0 <- est[3]
	a1 <- est[4]
	b1 <- est[5]
	c1 <- est[6]
	b0 <- a0 + (a1 - b1) * xj
	c0 <- b0 + (b1 - c1) * xk
	list(eta1 = eta1, eta2 = eta2, a0 = a0, a1 = a1, b0 = b0, b1 = b1, c0 = c0, c1
		 = c1)
}
