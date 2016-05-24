## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##' @title C implementation of the sinc function.
##' @param x argument
##' @return sinc(x)
sinc <- function(x) .Call(sinc_c, x)

##' @title C implementation of Zolotarev's function to the power 1-alpha.
##' @param x argument in [0,pi]
##' @param alpha parameter in (0,1]
##' @return sin(alpha*x)^alpha * sin((1-alpha)*x)^(1-alpha) / sin(x)
##' @author Martin Maechler
## NB: Is exported, documented, demo'ed, but not used by our own code!
A..Z <- function(x, alpha, I.alpha = 1 - alpha)
    .Call(A__c, x, alpha, I.alpha)

##' @title tan(pi*x), exact for integer x
##' @param x numeric vector
##' @return numeric vector of values tan(pi*x)
##' @author Martin Maechler
tanpi <- function(x) tan(pi * (x %% 1))

##' @title cos(pi/2 * x), exact for integer x
##' @param x numeric vector
##' @return numeric vector of values cos(pi/2 *x)
##' @author Martin Maechler
cospi2 <- function(x) {
    x <- r <- (x %% 4)## cos(pi/2 x) == cos(pi/2(x + 4k))  \forall k \in \Z
    if(any(isI <- x == round(x))) {
	i <- which(isI)
	r[i	      ] <-  0 # for those where x is 1 or 3
	r[i[x[i] == 0]] <-  1
	r[i[x[i] == 2]] <- -1
    }
    io <- which(!isI)
    r[io] <- cos(pi/2 * x[io])
    r
}

##' Sample S ~ S(alpha, beta, gamma, delta; pm), see package
##' \pkg{stabledist} for the parameterization.
##'
##' @title Sampling stable distributions
##' @param n number of random variates to be generated
##' @param alpha, see code in fBasics
##' @param beta, see code in fBasics
##' @param gamma, see code in fBasics
##' @param delta, see code in fBasics
##' @param pm in {0,1} parameterization, see code in fBasics
##' @return vector of variates S
##' @author Martin Maechler, based on Diethelm Wuertz's code in fBasics
rstable1R <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 1)
{
    stopifnot((la <- length(alpha)) >= 1, (lb <- length(beta)) >= 1,
	      length(gamma) >= 1, length(delta) >= 1,
	      0 < alpha, alpha <= 2, abs(beta) <= 1,
	      length(pm) == 1, pm %in% 0:1)

    p2 <- pi/2
    ## Special case (a,b) = (1,0):
    if (all(alpha == 1) && all(beta == 0)) {
	Z <- rcauchy(n)
    }
    else {
	## MM: Nolan(2009) "chapt1.pdf", p.21 has  "Theorem 1.19"
	## -- and attributes that to  'Chambers et al. (1976)'

	## Calculate uniform and exponential distributed random numbers:
	Theta <- pi * (runif(n)-1/2)
	W <- rexp(n)
	##  ^^^^^^ was "-log(runif(n))"
	## rexp() is faster, giving different numbers

        a.is.vec <- (la > 1)          # alpha is "vector" (not scalar)
        if(a.is.vec) {
            ## if alpha is not scalar, make sure that lengths of
            ##		alpha, beta, Theta, W  are all equal (== n)
            alpha <- rep(alpha, length.out = n)
            beta  <- rep(beta,  length.out = n)
        }

	norm <- alpha != 1 ## TODO:  abs(alpha - 1) > eps.alpha1
	## FIXME(2): ditto for	  | alpha - 1 | << 1

	Z <- numeric(n)
	if(any(norm)) { ## alpha != 1
	    alp <- alpha[norm]; Thet <- Theta[norm]
	    b.tan.pa <- beta[norm]*tanpi(alp/2)
	    th0 <- atan(b.tan.pa) / alp ## == \theta_0
	    ## Now, from Nolan/Chambers' formula, we replace
	    ##	1/(\cos(\alpha\theta_0) \cos\Theta)^{1/\alpha} with
	    ## c / (\cos\Theta)^{1/\alpha} where
	    ## c := (1 + (\beta*\tan(\pi\alpha/2))^2)^{1/{2\alpha}}
	    ## need to show that c = 1/(\cos(\alpha\theta_0))^{1/\alpha}
	    ## <==> 1 + (\beta*\tan(\pi\alpha/2))^2 = 1
            ## / (\cos(\alpha\theta_0))^2 and that's true, as
            ## 1 + (tan(y))^2 = 1 / (cos(y))^2 for
            ##   y = \alpha\theta_0 = \arc\tan(\beta*\tan(\pi\alpha/2))
	    c. <- (1 + b.tan.pa^2)^(1/(2*alp))
	    a.tht <- alp*(Thet+th0)
	    Z[norm] <-
		sin(a.tht) * c. / cos(Thet)^(1/alp) *
		    (cos(a.tht-Thet)/W[norm])^((1-alp)/alp)
	}
        ## {note that logicals vectorize, indices do *not* so easily}
	if(any(a1 <- !norm)) { ## alpha == 1
	    bet <- beta[a1]; Thet <- Theta[a1]
	    p2.bt <- p2 + bet*Thet
	    Z[a1] <- (p2.bt*tan(Thet) - bet*log((p2*W[a1]*cos(Thet))/p2.bt))/p2
	}
    }

    if(pm == 0)
        ## delta_1 := delta_0 - .. [Nolan, chapt.1, (1.7)],
        ## since above 1-parametr.
	delta <- delta - gamma * {
	    if(a.is.vec) {
		d.of <- numeric(n)
		d.of[norm] <- b.tan.pa
		d.of[a1	 ] <- beta * log(gamma)/p2
		d.of
	    }
	    else { ## alpha is scalar
		if(norm) b.tan.pa else beta * log(gamma)/p2
	    }
	}

    ## Result: Location - Scale trafo -- only now using (gamma, delta):
    Z * gamma + delta
}

##' Sample S ~ S(alpha, beta, gamma, delta; pm), see rstable1R() above.
##' For beta == 1 and pm == 1, the fast C implementation is used.
##'
##' @title Efficiently sampling stable distributions
##' @param n number of random variates to be generated
##' @param alpha, see code in fBasics
##' @param beta, see code in fBasics
##' @param gamma, see code in fBasics
##' @param delta, see code in fBasics
##' @param pm in {0,1} parameterization, see code in fBasics
##' @return vector of variates S
##' @author Martin Maechler
rstable1C <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 1)
{
    stopifnot((la <- length(alpha)) >= 1, (lb <- length(beta)) >= 1,
	      length(gamma) >= 1, length(delta) >= 1,
	      0 < alpha, alpha <= 2, abs(beta) <= 1, gamma >= 0,
              length(pm) == 1, pm %in% 0:1)

    if(beta == 1 && pm == 1)
        .Call(rstable_c, n, alpha) * gamma + delta
    else rstable1R(n, alpha=alpha, beta=beta,
                   gamma=gamma, delta=delta, pm=pm)
}

## Now works (but gains almost no speed !)
rstable1 <- rstable1C
