negbin1 <- function (link, phi = stop("'phi' must be given"))
{
	.Phi <- phi
	env <- new.env(parent = .GlobalEnv)
	assign(".Phi", phi, envir = env)
	stats <- make.link("identity")
	variance <- function(mu) mu * (1 + .Phi)
  valideta <- function(eta) all(eta >= 0)
	validmu <- function(mu) all(mu >= 0)
	dev.resids <- function(y, mu, wt) {
		r <- (mu - y)/.Phi * log(1 + .Phi)
		p <- which(y > 0)
		r[p] <- r[p] + (mapply(function(Y,m,phi) {
								k <- seq_len(Y) - 1
								sum(log((Y/phi + k) / (m/phi + k)))},
							y[p], mu[p], MoreArgs = list(phi = .Phi)))
		2 * r
	}
	aic <- function(y, n, mu, wt, dev) -2 * sum(dnbinom(y, size = mu/.Phi, prob = 1/(1+.Phi), log = TRUE))
	initialize <- expression({
		if (any(y < 0)) stop("negative values not allowed for the 'negbin1' family")
		if (any(abs(y - round(y)) > 0.001)) stop("non-integer counts in a negbin model")
		n <- rep.int(1, nobs)
	})
	simfun <- function(object, nsim) {
		ftd <- fitted(object)
		val <- rnbinom(nsim * length(ftd), size = ftd/.Phi, prob = 1/(1+.Phi))
		val
	}
    fname <- paste("negbin1 (phi = ", format(round(phi, 4)), ")", sep = "")
	environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
	structure(list(family = fname, link = "identity", linkfun = stats$linkfun,
		linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
		aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
		valid.mu = validmu, valid.eta = valideta, simulate = simfun),
		class = "family")
}