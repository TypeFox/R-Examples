##  TWEEDIEF.R

tweedie <- function(var.power=0, link.power=1-var.power) {
#	Tweedie generalized linear model family
#	Gordon Smyth
#	22 Oct 2002.  Last modified 2 Sep 2011.

	lambda <- link.power
	if(lambda==0) {
		linkfun <- function(mu) log(mu)
		linkinv <- function(eta) pmax(.Machine$double.eps, exp(eta))
		mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
		valideta <- function(eta) TRUE
	} else {
		linkfun <- function(mu) mu^lambda
		linkinv <- function(eta) eta^(1/lambda)
		mu.eta <- function(eta) (1/lambda) * eta^(1/lambda - 1)
		valideta <- function(eta) TRUE
	}
	p <- var.power
	variance <- function(mu) mu^p
	if(p == 0)
		validmu <- function(mu) TRUE
	else if(p > 0)
		validmu <- function(mu) all(mu >= 0)
	else
		validmu <- function(mu) all(mu > 0)
	dev.resids <- function(y, mu, wt) {
		y1 <- y + 0.1*(y == 0)
		if (p == 1)
			theta <- log(y1/mu)
		else
			theta <- ( y1^(1-p) - mu^(1-p) ) / (1-p)
		if (p == 2)
#			Returns a finite somewhat arbitrary residual for y==0, although theoretical value is -Inf
			kappa <- log(y1/mu)
		else
			kappa <- ( y^(2-p) - mu^(2-p) ) / (2-p)
		2 * wt * (y*theta - kappa)
	}	
	initialize <- expression({
		n <- rep(1, nobs)
		mustart <- y + 0.1 * (y == 0)
	})
	aic <- function(y, n, mu, wt, dev) NA
	structure(list(
		family = "Tweedie", variance = variance, dev.resids = dev.resids, aic = aic,
		link = paste("mu^",as.character(lambda),sep=""), linkfun = linkfun, linkinv = linkinv,
		mu.eta = mu.eta, initialize = initialize, validmu = validmu, valideta = valideta), 
		class = "family")
}

