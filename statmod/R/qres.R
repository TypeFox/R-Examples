## QRES.R

qresiduals <- qresid <- function(glm.obj, dispersion=NULL)
#	Wrapper function for quantile residuals
#	Peter K Dunn
#	28 Sep 2004.  Last modified 5 Oct 2004.
{
glm.family <- glm.obj$family$family
if(substr(glm.family,1,17)=="Negative Binomial") glm.family <- "nbinom"
switch(glm.family,
   binomial = qres.binom( glm.obj),
   poisson = qres.pois(glm.obj),
   Gamma = qres.gamma(glm.obj, dispersion),
   inverse.gaussian = qres.invgauss(glm.obj, dispersion),
   Tweedie = qres.tweedie(glm.obj, dispersion),
   nbinom = qres.nbinom(glm.obj),
   qres.default(glm.obj, dispersion)
)}

qres.binom <- function(glm.obj)
#	Randomized quantile residuals for binomial glm
#	Gordon Smyth
#	20 Oct 96.  Last modified 25 Jan 02.
{
	p <- fitted(glm.obj)
	y <- glm.obj$y
	if(!is.null(glm.obj$prior.weights))
		n <- glm.obj$prior.weights
	else
		n <- rep(1,length(y))
	y <- n * y
	a <- pbinom(y - 1, n, p)
	b <- pbinom(y, n, p)
	u <- runif(n = length(y), min = a, max = b)
	qnorm(u)
}

qres.pois <- function(glm.obj)
#	Quantile residuals for Poisson glm
#	Gordon Smyth
#	28 Dec 96
{
	y <- glm.obj$y
	mu <- fitted(glm.obj)
	a <- ppois(y - 1, mu)
	b <- ppois(y, mu)
	u <- runif(n = length(y), min = a, max = b)
	qnorm(u)
}

qres.gamma <- function(glm.obj, dispersion = NULL)
#	Quantile residuals for gamma glm
#	Gordon Smyth
#	28 Dec 96.  Last modified 10 Jan 97
{
	mu <- fitted(glm.obj)
	y <- glm.obj$y
	df <- glm.obj$df.residual
	w <- glm.obj$prior.weights
	if(is.null(w)) w <- 1
	if(is.null(dispersion)) dispersion <- sum(w * ((y - mu)/mu)^2)/df
	u <- pgamma((w * y)/mu/dispersion, w/dispersion)
	qnorm(u)
}

qres.invgauss <- function(glm.obj, dispersion = NULL)
#	Quantile residuals for inverse Gaussian glm
#	Gordon Smyth
#	Created 15 Jan 98. Last modified 31 May 2014.
{
	mu <- fitted(glm.obj)
	y <- glm.obj$y
	df <- glm.obj$df.residual
	w <- glm.obj$prior.weights
	if(is.null(w)) w <- 1
	if(is.null(dispersion)) dispersion <- sum(w * (y - mu)^2 / (mu^2*y)) / df
	up <- y>mu
	down <- y<mu
	if(any(down)) y[down] <- qnorm(pinvgauss(y,mean=mu,dispersion=dispersion,log.p=TRUE),log.p=TRUE)
	if(any(up)) y[up] <- qnorm(pinvgauss(y,mean=mu,dispersion=dispersion,lower.tail=FALSE,log.p=TRUE),lower.tail=FALSE,log.p=TRUE)
	y
}

qres.nbinom <- function(glm.obj)
{
#	Quantile residuals for Negative Binomial glm
#	Gordon Smyth
#	22 Jun 97.  Last modified 18 Nov 2008.
#
	y <- glm.obj$y
	if(is.null(glm.obj$theta)) {
		size <- glm.obj$call$family[[2]]
	} else {
		size <- glm.obj$theta
	}
	mu <- fitted(glm.obj)
	p <- size/(mu + size)
	a <- ifelse(y > 0, pbeta(p, size, pmax(y, 1)), 0)
	b <- pbeta(p, size, y + 1)
	u <- runif(n = length(y), min = a, max = b)
	qnorm(u)
}

qres.tweedie <- function(glm.obj, dispersion = NULL)
#	Quantile residuals for Tweedie glms
#	Gordon Smyth
#	Created 29 April 1998.  Last modified 30 March 2015.
{
	requireNamespace("tweedie")
	mu <- fitted(glm.obj)
	y <- glm.obj$y
	df <- glm.obj$df.residual
	w <- glm.obj$prior.weights
	if(is.null(w))
		w <- 1
	p <- get("p",envir=environment(glm.obj$family$variance))
	if(is.null(dispersion))
		dispersion <- sum((w * (y - mu)^2)/mu^p)/df
	u <- tweedie::ptweedie(q=y, power=p, mu=fitted(glm.obj), phi=dispersion/w)
	if(p>1&&p<2)
		u[y == 0] <- runif(sum(y == 0), min = 0, max = u[y == 0])
	qnorm(u)
}

qres.default <- function(glm.obj, dispersion=NULL)
#	Quantile residuals for Gaussian and default glms
#	Gordon Smyth
#	5 Oct 2004.
{
	r <- residuals(glm.obj, type="deviance")
	if(is.null(dispersion)) {
		df.r <- glm.obj$df.residual
		if(df.r > 0) {
			if(any(glm.obj$weights==0)) warning("observations with zero weight ", "not used for calculating dispersion")
	        dispersion <- sum(glm.obj$weights*glm.obj$residuals^2)/df.r
	    } else
	    	dispersion <- 1
    }
    r/sqrt(dispersion)
}

