Digamma <- function(link = "log") {
#	Digamma generalized linear model family
#	Gordon Smyth, smyth@wehi.edu.au
#  3 July 1998.  Last revised 9 Dec 2002.
#
#	improve on the link deparsing code in quasi()
	linkarg <- substitute(link)
	if (is.expression(linkarg) || is.call(linkarg)) {
		linkname <- deparse(linkarg)
	} else if(is.character(linkarg)) { 
		linkname <- linkarg
		link <- make.link(linkarg)
	} else if(is.numeric(linkarg)) {
		linkname <- paste("power(",linkarg,")",sep="")
		link <- make.link(linkarg)
	} else {
		linkname <- deparse(linkarg) 
		link <- make.link(linkname)
	}
	validmu <- function(mu) all(mu>0)
	dev.resids <- function(y, mu, wt) wt * unitdeviance.digamma(y,mu)
	initialize <- expression({
		if (any(y <= 0)) stop(paste("Non-positive values not", "allowed for the Digamma family"))
		n <- rep(1, nobs)
		mustart <- y
	})
	aic <- function(y, n, mu, wt, dev) NA
	structure(list(
		family = "Digamma",
		variance = varfun.digamma, dev.resids = dev.resids, aic = aic,
		link = linkname,
		linkfun = link$linkfun, linkinv = link$linkinv, mu.eta = link$mu.eta,
		valideta = link$valideta, validmu = validmu, initialize = initialize, 
		class = "family"))
}

cumulant.digamma <- function(theta)
#	Cumulant function for the Digamma family
#	GKS  3 July 98
	2*( theta*(log(-theta)-1) + lgamma(-theta) )

meanval.digamma <- function(theta)
#	Mean value function for the Digamma family
#	GKS  3 July 98
	2*( log(-theta) - digamma(-theta) )

d2cumulant.digamma <- function(theta)
#	2nd derivative of cumulant function for Digamma family
#	GKS  3 July 98
	2*( 1/theta + trigamma(-theta) )

canonic.digamma <- function(mu) {
#	Canonical mapping for Digamma family
#	Solve meanval.digamma(theta) = mu for theta
#	GKS  3 July 98
#
#	Starting value from -log(-theta) =~ log(mu)
	mlmt <- log(mu)	
	theta <- -exp(-mlmt)

	for (i in 1:3) {
		mu1 <- meanval.digamma(theta)
		v <- d2cumulant.digamma(theta)
		deriv <- -v/mu1*theta
		mlmt <- mlmt - log(mu1/mu)/deriv	
		theta <- -exp(-mlmt)
	}
	theta
}

varfun.digamma <- function(mu) {
#	Variance function for Digamma family
#	GKS  3 July 98
#
	theta <- canonic.digamma(mu)
	2*( 1/theta + trigamma(-theta) )
}

unitdeviance.digamma <- function(y,mu) {
#	Unit deviance for Digamma family
#	GKS  3 July 98
#
	thetay <- canonic.digamma(y)
	theta <- canonic.digamma(mu)
	2*( y*(thetay-theta) - (cumulant.digamma(thetay)-cumulant.digamma(theta)) )
}
