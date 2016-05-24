# estimate the Bayesian posterior probability of each alternative being the best poisson bandit
best_poisson_bandit <-
function(x, n=NULL) {
	k <- length(x)
  # gamma, using 'non-informative' prior alpha=1, beta=0
  # posterior is alpha = 1 + sum(x[[i]]), beta = length(x[[i]])
	alphas = x
	betas = n
	if (is.null(n)) {
  	alphas = sapply(x, sum)
  	betas = sapply(x, length)	  
	}
	means = alphas/betas
	alphas = alphas + 1

  # just integrate up to a max range -- because integrate doesn't work well if the function is near 0 over most of its support
  max_range = max(means) * 2

	ans <- numeric(k)
	for (i in (1:k)) {
		indx <- (1:k)[-i]
		f <- function(z) {
		  # probability that this outcome's parameter is at this point and all other outcomes' parameters are lower
			r <- dgamma(z, shape=alphas[i], rate=betas[i])
			for (j in indx) {
				r <- r * pgamma(z, shape=alphas[j], rate=betas[j])
			}
			return(r)
		}
		ans[i] = integrate(f,lower=0,upper=max_range,subdivisions=10000)$value
	}
	return(ans)
}

bpb <-
function(x) {
	best_poisson_bandit(x)
}
