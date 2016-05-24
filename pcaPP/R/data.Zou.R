
data.Zou <- function (n = 250, p =  c(4, 4, 2), ...)
{	##	generating testdata as Zou 2004

	r.var <- c (290, 300, 1)
	V <- matrix (rnorm (n * 3), ncol = 3) %*% diag (sqrt (r.var))

	V[,3] <- -0.3 * V[,1] + 0.925 * V[,2] + V[,3]

	dat <- matrix (c(rep (V[,1], p[1]), rep (V[,2], p[2]), rep (V[,3], p[3])), nrow = n)
	err <- rnorm (length (dat))

	dat + err
}
