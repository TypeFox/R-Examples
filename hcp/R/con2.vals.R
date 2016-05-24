con2.vals <-
function(x, y, n, j, k)
{
	a <- con2.parms(x, y, n, j, k, 1, 1)
	nr <- nrow(a$theta)
	est <- a$theta[nr,  ]
	b <- con2.est(x[j], x[k], est)
	eta <- c(b$eta1, b$eta2)
	beta <- c(b$a0, b$a1, b$b0, b$b1, b$b2, b$c0, b$c1)
	tau <- c(x[j], x[k])
	list(eta = eta, beta = beta, tau = tau)
}
