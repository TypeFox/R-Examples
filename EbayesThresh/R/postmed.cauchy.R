"postmed.cauchy" <-
function(x, w)
{
#
# find the posterior median of the Cauchy prior with
#   mixing weight w, pointwise for each of the data points x
#
	nx <- length(x)
	zest <- rep(NA, length(x))
	w <- rep(w, length.out = nx)
	ax <- abs(x)
	j <- (ax < 20)
	zest[!j] <- ax[!j] - 2/ax[!j]
	if(sum(j) > 0) {
		zest[j] <- vecbinsolv(zf = rep(0, sum(j)), fun = cauchy.medzero,
			tlo = 0, thi = max(ax[j]), z = ax[j], w = w[j])
	}
	zest[zest < 1e-007] <- 0
	zest <- sign(x) * zest
	return(zest)
}
