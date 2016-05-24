"tccfAR" <-
function(phi, theta)
{
#auxilary function used with iarma#########
#computes the theoretical cross-covariance function of two autoregressions
# z[t]-phi[1] z_[t-1] --- phi[p] z[t-p]     = a[t]
# z[t]-theta[1] z_[t-1] --- theta[q] z[t-q] = a[t]
# where p, q are length(phi), length(theta)
	p <- length(phi)
	q <- length(theta)
	if(p == 0 || q == 0)
		return(numeric(0))
	k <- p + q
	rhs <- c(-1, rep(0, k - 1))
	A <- matrix(numeric(k^2), nrow = k, ncol = k)
	for(i in 1:k) {
		for(j in 1:k) {
			imj <- i - j
			ijq <- i + j - q - 1
			if(i > q) {
				if(i > j && imj <= q)
				  A[i, j] <- theta[imj]
				else if(i > q && imj == 0)
				  A[i, j] <- -1
			}
			else {
				if(ijq > 0 && ijq <= p)
				  A[i, j] <- phi[ijq]
				else if(ijq == 0)
				  A[i, j] <- -1
			}
		}
	}
	return(solve(A, rhs))
}

