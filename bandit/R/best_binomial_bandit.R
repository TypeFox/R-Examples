# estimate the Bayesian posterior probability of each alternative being the best binomial bandit
best_binomial_bandit <-
function(x, n, alpha=1, beta=1) {
	k <- length(x)
	ans <- numeric(k)
	for (i in (1:k)) {
			indx <- (1:k)[-i]
			f <- function(z) {
					r <- dbeta(z, x[i] + alpha, n[i] - x[i] + beta)
					for (j in indx) {
							r <- r * pbeta(z, x[j] + alpha, n[j] - x[j] + beta)
					}
					return(r)
			}
			ans[i] = integrate(f, 0, 1)$value
	}
	return(ans)
}

bbb <-
function(x, n, alpha=1, beta=1) {
	best_binomial_bandit(x, n, alpha, beta)
}
