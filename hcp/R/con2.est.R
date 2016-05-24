con2.est <-
function(xj, xk, est)
{
	eta1 <- est[1]
	eta2 <- est[2]
	a0 <- est[3]
	a1 <- est[4]
	b1 <- est[5]
        b2 <- est[6]
	c1 <- est[7]
	b0 <- a0 + (a1 - b1) * xj - b2 * xj^2
	c0 <- b0 + (b1 - c1) * xk + b2 * xk^2
	list(eta1 = eta1, eta2 = eta2, a0 = a0, a1 = a1, b0 = b0, b1 = b1, 
             b2 = b2, c0 = c0, c1 = c1)
}
