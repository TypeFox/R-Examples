### Lévy :
### ====
## Stable(alpha = 1/2, beta = 1, gamma, delta, pm = 1)	<===>  Levy(delta, gamma)
##	  ~~~~~~~~~~~  ~~~~~~~~
## http://en.wikipedia.org/wiki/L%C3%A9vy_distribution
## The probability density function of the Lévy distribution over the domain x >= \mu is
##
##     f(x;\mu,c)=\sqrt{\frac{c}{2\pi}}~~\frac{e^{ -\frac{c}{2(x-\mu)}}} {(x-\mu)^{3/2}}
##
## NOTA BENE: You can use 'mpfr numbers for x' -- ! --> ../../tests/pstab-ex.R
##							~~~~~~~~~~~~~~~~~~~~~~
dLevy <- function(x, mu=0, c=1, log=FALSE) {
    r <- x <- x-mu
    ## ensure f(0) = 0 {not NaN}:
    pos <- x > 0 ; x <- x[pos]; if(log) r[!pos] <- -Inf
    r[pos] <- if(log)
	(log(c/(2*pi)) + -c/x - 3*log(x))/2
    else
	sqrt(c/(2*pi)) * exp(-c/(2*x)) / (x^(3/2))
    r
}
## where \mu is the location parameter and c is the scale parameter.

## The cumulative distribution function is
##
##     F(x;\mu,c)=\textrm{erfc}\left(\sqrt{c/(2(x-\mu))}\right)
##
##     {MM:  fixed Wikipedia entry: (x-mu) is in the denominator!}
pLevy <- function(x, mu=0, c=1, log.p=FALSE, lower.tail=TRUE) {
    ## erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
    ## erfc(sqrt(c/(2*(x-mu))))
    x <- (x-mu)/c # re-scale to (0,1)
    u <- 1/sqrt(x)
    if(log.p) {
	if(lower.tail)
	    log(2) + pnorm(u, lower.tail = FALSE, log.p=TRUE)
	else log(2 * pnorm(u) - 1)
    } else {
	if(lower.tail)
	    2* pnorm(u, lower.tail = FALSE)
	else 2*pnorm(u) - 1
    }
}

## where \textrm{erfc}(z) is the complementary error function. The shift
## parameter \mu has the effect of shifting the curve to the right by an
## amount \mu, and changing the support to the interval [\mu, \infty). Like
## all stable distributions, the Levy distribution has a standard form
## f(x;0,1) which has the following property:
##
##     f(x;\mu,c) dx  =	 f(y;0,1) dy
##
## where y is defined as
##
##     y = \frac{x-\mu}{c}
