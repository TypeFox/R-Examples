varbin <-
function(n, m, alpha = 0.05, R = 5000){
	
	dat <- data.frame(n = n, m = m)

  n <- dat$n; m <- dat$m; N <- length(n); ntot <- sum(n); mtot <- sum(m)
  
	# binomial
  mu <- mtot / ntot
  varmu <- mu * (1 - mu) / (ntot - 1)
  
	# ratio
  muratio <- mu
  varmuratio <- N * (N - 1)^(-1) * ntot^(-2) * sum((m - n * mu)^2)
  
	# arithmetic
  muarithm <- mean(m / n)
  varmuarithm <- var(m / n) / N
  
	# jackknife
  z <- (m - n * mu) / (ntot - n)
  pseudovalue <- mu + (N - 1) * z
  #mujack <- mu + (N - 1) * mean(z)
  #varmujack <- (N^(-1)) * (N - 1) * sum((z - mean(z))^2)
  mujack <- mean(pseudovalue)
  varmujack <- var(pseudovalue) / N
  
  # bootstrap
  if(!require(boot, quietly = TRUE))
    stop("This function requires the recommended package dQuote(boot).")
  foo <- function(d, f) mu <- sum(d$m * f) / sum(d$n * f)
  res <- boot(data = dat, statistic = foo, stype = "w", R = R)
  muboot <- mean(res$t)
  varmuboot <- var(res$t)
  
	# results
  tab <- data.frame(mu = c(mu, muratio, muarithm, mujack, muboot),
                    varmu = c(varmu, varmuratio, varmuarithm, varmujack, varmuboot))
  tab$lower <- pmax(0, tab$mu - qnorm(1 - alpha/2) * sqrt(tab$varmu))
  tab$upper <- pmin(1, tab$mu + qnorm(1 - alpha/2) * sqrt(tab$varmu))
  rownames(tab) <- c("Binomial", "Ratio", "Arithmetic", "Jackknife", "Bootstrap")
  features <- c(N = N, n = ntot, m = mtot)
  
	# outputs
  structure(
  	list(tab = tab, muboot = res$t, alpha = alpha, features = features, dat = dat),
  	class = "varbin"
  	)

}
