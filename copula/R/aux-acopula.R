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


if(getRversion() < "2.15")
paste0 <- function(...) paste(..., sep="")

#### Functions and Methods for "acopula" objects
#### class definition in ./AllClass.R

### Ali-Mikhail-Haq == "AMH" ###################################################

##' @title Ali-Mikhail-Haq ("AMH")'s  tau(theta)
##' @param theta
##' @return 1 - 2*((1-th)*(1-th)*log(1-th)+th)/(3*th*th)  where th := theta
##' numerically accurately, for both limits  th -> 0  &  th -> 1
##' Nelsen (2006, p.172) range of tau: [(5 - 8 log 2) / 3, 1/3] ~= [-0.1817, 0.3333]
##' @author Martin Maechler
tauAMH <- function(theta) {
    if(length(theta) == 0) return(numeric(0)) # to work with NULL
    r <- theta
    na <- is.na(theta)
    lrg <- (abs(theta) > 0.01) & !na ## for "large" theta (<=> theta >= 0.01), use standard formula
    f. <- function(t) {
	## 1 - 2*((1-t)*(1-t)*log1p(-t) + t) / (3*t*t)
	r <- t
	r[i1 <- (1-t) == 0] <- 1/3
        t <- t[s <- !i1]
	r[s] <- 1 - 2*((1-t)^2 *log1p(-t) + t) / (3*t*t)
	r
    }
    r[lrg] <- f.(theta[lrg])
    ## Otherwise use "smart" Taylor polynomials: for given order, *drop* the last two..
    ## with the following coefficients:
    ##_ x <- polynomial()
    ##_ MASS::fractions(as.vector(mkAMHtau.0series(12,TRUE) * 9/(2*x)))
    ##  1  1/4  1/10  1/20  1/35  1/56  1/84  1/120  1/165  1/220

    if(any(!lrg & !na)) {
	l1 <- !lrg & !na & (ll1 <- theta > 2e-4) ## --> k = 7
	r[l1] <- (function(x) 2*x/9*(1+ x*(1/4 + x/10*(1 + x*(1/2 + x/3.5))))
		  )(theta[l1])
	l2 <- !ll1 & !na & (ll2 <- theta > 1e-5)	 ## --> k = 6
	r[l2] <- (function(x) 2*x/9*(1+ x*(1/4 + x/10*(1 + x/2))))(theta[l2])
	l3 <- !ll2 & !na & (ll <- theta > 2e-8)	## --> k = 5
	r[l3] <- (function(x) 2*x/9*(1+ x*(1/4 + x/10)))(theta[l3])
	irem <- which(!ll & !na)## rest: theta <= 2e-8 : k == 4
	r[irem] <- (function(x) 2*x/9*(1+ x/4))(theta[irem])
    }
    r
}


### Clayton ####################################################################

##' Note: this is mainly to show that this function can be very well
##' approximated much more simply by just using m <- round(V0).
##'
##' @title Optimal constant for fast rejection
##' @param V0 numeric vector >= 0
##' @return optimal constant m for the fast rejection algorithm
##' @author Martin Maechler (based on Marius Hofert's code)
m.opt.retst <- function(V0) {
    n <- length(V0)
    fV <- floor(V0)
    cV <- ceiling(V0)
    v1 <- fV*exp(V0/fV)
    v2 <- cV*exp(V0/cV)

    m <- integer(n)
    l1 <- (V0 <= 1)
    m[which(l1)] <- 1L

    i2 <- which(!l1) ## those with V0 > 1
    l3 <- (v1[i2] <= v2[i2])
    i3 <- i2[l3]
    m[i3] <- fV[i3]
    i4 <- i2[!l3]
    m[i4] <- cV[i4]
    m
}

### fast rejection for fixed m, R version

##' Sample a random variate St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((1+t)^alpha-1)), see Nolan's book for
##' the parametrization, via an m-fold sum of random variates from
##' \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0/m)^{1/alpha}, (V_0/m)
##' *I_{alpha = 1}, I_{alpha != 1}; 1) with Laplace-Stieltjes transform
##' exp(-(V_0/m)*((1+t)^alpha-1)). This is a building block for the fast rejection.
##'
##' @title Sample an exponentially tilted stable distribution as an m-fold sum
##' @param m number of summands, any positive integer
##' @param V0 random variate
##' @param alpha parameter in (0,1]
##' @return St
##' @author Marius Hofert, Martin Maechler
retstablerej <- function(m,V0,alpha) {
    gamm. <- (cospi2(alpha)*V0/m)^(1/alpha)
    sum(unlist(lapply(integer(m),
		      function(.) {
			  ## apply standard rejection for sampling
			  ## \tilde{S}(alpha, 1, (cos(alpha*pi/2)
			  ##	*V_0/m)^{1/alpha}, (V_0/m)*I_{alpha = 1},
			  ## h*I_{alpha != 1}; 1) with Laplace-Stieltjes
			  ## transform exp(-(V_0/m)*((h+t)^alpha-h^alpha))
			  repeat {
			      V__ <- rstable1(1, alpha, beta=1, gamma = gamm.)
			      if(runif(1) <= exp(-V__))
				  return(V__)
			  }})
	       ## on acceptance, St_k ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
	       ## *V_0/m)^{1/alpha}, (V_0/m)*I_{alpha = 1}, h*I_{alpha != 1};
	       ## 1) with Laplace-Stieltjes transform
	       ## exp(-(V_0/m)*((h+t)^alpha-h^alpha))
	       ))
}

### fast rejection algorithm, R version

##' Sample a vector of random variates St ~ \tilde{S}(alpha, 1,
##' (cos(alpha*pi/2)*V_0)^{1/alpha}, V_0*I_{alpha = 1},
##' h*I_{alpha != 1}; 1) with LS transform
##' exp(-V_0((h+t)^alpha-h^alpha)) with the fast rejection
##' algorithm; see Nolan's book for the parametrization
##'
##' @title Sampling an exponentially tilted stable distribution
##' @param alpha parameter in (0,1]
##' @param V0 vector of random variates
##' @param h non-negative real number
##' @return vector of variates St
##' @author Marius Hofert, Martin Maechler
retstableR <- function(alpha, V0, h=1) {
    stopifnot(is.numeric(alpha), length(alpha) == 1,
	      0 <= alpha, alpha <= 1) # alpha > 1 => cos(pi/2 *alpha) < 0
    n <- length(V0)
    ## case alpha == 1
    if(alpha == 1 || n == 0) return(V0) # alpha == 1 => point mass at V0
    ## else alpha != 1 => call fast rejection algorithm with optimal m
    m <- m.opt.retst(V0)
    mapply(retstablerej, m=m, V0=V0, alpha=alpha)
}

### state-of-the-art: fast rejection + Luc's algorithm, C version

##' Sample random variates St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)), see Nolan's book for
##' the parametrization, with the fast rejection.
##' This procedure is more efficient than retstableR since it calls the C
##' function retstable_c and uses both the fast rejection and Luc Devroye's algorithm.
##'
##' @title Efficiently sampling an exponentially tilted stable distribution
##' @param alpha parameter in (0,1]
##' @param V0 vector of random variates
##' @param h non-negative real number
##' @param method which method to call ("Marius Hofert", "Luc Devroye")
##' @return vector of variates St
##' @author Martin Maechler
retstableC <- function(alpha, V0, h = 1, method = NULL) {
    n <- length(V0)
    stopifnot(is.numeric(alpha), length(alpha) == 1,
	      0 < alpha, alpha <= 1,
	      is.numeric(h), length(h) == 1, h > 0)
    if(alpha == 1 || n == 0) {
	## alpha == 1 => St corresponds to a point mass at V0 with
	V0           # Laplace-Stieltjes transform exp(-V0*t)
    }
    else {
	if(is.null(method)) {
	    if(any(diff(V.is.sml <- V0 * h^alpha < 4))) { ## use *both* methods
		r <- numeric(n)
		r[ V.is.sml] <- .Call(retstable_c, V0[ V.is.sml], h = h, alpha, "MH")
		r[!V.is.sml] <- .Call(retstable_c, V0[!V.is.sml], h = h, alpha, "LD")
		return(r)
	    }
	    else
		method <- if(V.is.sml[1]) "MH" else "LD"
	}
	else
	    method <- match.arg(method, c("MH","LD"))
	.Call(retstable_c, V0, h = h, alpha, method)
    }
}

## switch to make fast C code the default
if(FALSE)
    retstable <- retstableR
retstable <- retstableC # retstable is by default retstableC


### Frank ######################################################################

### sampling a logarithmic distribution, R version

##' Random number generator for a Log(p) distribution with the algorithm "LK" of
##' Kemp (1981), R version.
##'
##' @title Sample a Log(p) distribution
##' @param n number of random variates to be generated
##' @param p parameter in (0,1)
##' @param Ip = 1 - p_ (possibly more accurate) -- use, if p ~= 1
##' @return vector of random variates from Log(p)
##' @author Marius Hofert, Martin Maechler
rlogR <- function(n, p, Ip = 1-p) {
    if(missing(p)) p <- 1 - Ip
    stopifnot((n <- as.integer(n)) >= 0,
              0 < p, p <= 1, 0 < Ip)
    vec <- numeric(n)
    if(n >= 1) {
	u <- runif(n)
	l1 <- u > p
	vec[l1] <- 1
	i2 <- which( !l1 ) # of shorter length, say n2
	q2 <- 1-(Iq2 <- Ip^runif(length(i2))) # length n2
	l3 <- u[i2] < q2*q2
	i3 <- i2[l3]
	vec[i3] <- floor(1+abs(log(u[i3])/log1p(-Iq2[l3])))
	l4 <- u[i2] > q2
	vec[i2[l4]] <- 1
	l5 <- ! (l3 | l4) # q2^2 <= u[i2] <= q2
	vec[i2[l5]] <- 2
    }
    vec
}

### state-of-the art: sampling a logarithmic distribution, C version

##' Random number generator for a Log(p) distribution with the algorithm "LK" of
##' Kemp (1981), C version.
##'
##' @title Efficiently sampling a Log(p) distribution
##' @param n number of random variates to be generated
##' @param p parameter in (0,1)
##' @param Ip = 1 - p_ (possibly more accurate)
##' @return vector of random variates from Log(p)
##' @author Martin Maechler
rlog <- function(n, p, Ip = 1-p) {
    if(missing(p)) p <- 1 - Ip
    stopifnot(n >= 0, 0 < p, p <= 1, 0 < Ip)
    .Call(rLog_vec_c, n, p, Ip)
}

### state-of-the art: sampling F01Frank, C version

##' Generate a vector of variates V01 ~ F01 with Laplace-Stieltjes transform
##' ((1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)))^V0.
##'
##' @title Efficiently sampling V01 for Frank
##' @param V0 vector of random variates from F0
##' @param theta0 parameter theta0 in (0,infinity)
##' @param theta1 parameter theta1 in [theta0, infinity)
##' @param rej method switch: if V0*theta_0*p0^(V0-1) > rej a rejection
##'        from F01 of Joe is applied (otherwise, the sum is
##'        sampled via a logarithmic envelope for the summands)
##' @param approx largest number of summands before asymptotics is used
##' @return vector of random variates V01
##' @author Marius Hofert
rF01Frank <- function(V0, theta0, theta1, rej, approx) {
    .Call(rF01Frank_vec_c, V0, theta0, theta1, rej, approx)
}

### wrapper for inner distribution F for Frank

##' Generate a vector of variates V ~ F with Laplace-Stieltjes transform
##' (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)).
##'
##' @title Sampling F for Frank
##' @param n number of variates from F
##' @param theta0 parameter theta0 in (0,infinity)
##' @param theta1 parameter theta1 in [theta0, infinity)
##' @param rej method switch: if theta_0 > rej a rejection
##'        from Joe (Sibuya distribution) is applied (otherwise, a logarithmic
##' 	   envelope is used)
##' @return vector of random variates V
##' @author Marius Hofert
rFFrank <- function(n, theta0, theta1, rej)
    rF01Frank(rep(1,n),theta0,theta1,rej,1) # approx = 1 (dummy)

## see documentation for more information
dDiagFrank <- function(u, theta, d, log=FALSE,
		       method = c("auto", paste0("poly",1:4),
                       "polyFull",  "m1", "MH", "MMH"))
{
    stopifnot(length(d <- as.integer(d)) == 1, d >= 1)
    r <- ut <- u*theta
    Ie <- -expm1(-ut) # 1 - exp(-u * theta)
    ## L <- ut > log(2)*.Machine$double.digits # <==> exp(-ut) < Mach..eps
    method <- match.arg(method)
    if(method == "auto")
	method <- {
	    if(d <= 3) "poly1" else # d in {1,2,3}
	    if(d <= 6) paste0("poly", d-2) # d in {4,5,6}
	    else ## this is not really good: but (u*th) is *vector*
		"polyFull"
    }
    if(substr(method, 1,4) == "poly") {
	e.ut <- exp(-ut)
	ep <- (e.ut - exp(ut-theta))/Ie # "FIXME": maybe improve, testing  u > 1/2
	d1 <- d-1
	D <- d + d1*ep # D = d + (d-1) * (exp(-theta*u)-exp(theta*u-theta))
        if(d > 2) { # <==> d1 = d-1 >= 2
            f <- 1 + ep # 1+exp(-theta*u)-exp(theta*u-theta)
            delt <- e.ut * f # exp(-theta*u) * (1+exp(-theta*u)-exp(theta*u-theta))
	    D <- D + f* delt *
		switch(method,
		       "poly1" = d1*(d-2L)/2,
		       "poly2" = d1*(d-2L)/2 *(1 + (d-3L)/3 * delt),
		       "poly3" = d1*(d-2L)/2 *(1 + (d-3L)/3 * delt *
                       (1 + (d-4L)/4 * delt)),
		       "poly4" = d1*(d-2L)/2 *(1 + (d-3L)/3 * delt *
                       (1 + (d-4L)/4 * delt*(1 + (d-5L)/5 * delt))),
		       "polyFull" = { ## full polynomial formula
			   if(is(ut, "mpfr"))
			       sfsmisc::polyn.eval(Rmpfr::chooseMpfr.all(d1)[-1],
                                                   x = delt)
			   else		 polynEval(choose(d1, 2:d1), x = delt)
		       },
		       stop("invalid poly* method: ", method,
			    ". Should never happen") )
	}
        if(log) log(d) - log(D) else d / D
    }
    else switch(method,
	   "m1" = {
	       h <- -expm1(-theta) # =	1 - exp(-theta)
	       D <- (h/Ie)^(d-1) - Ie # D := [d]enominator
	       if(log) log(d) -ut -log(D) else d*exp(-ut) / D
	   },
	   "MH" = { # Marius Hofert
	       h <- -expm1(-theta) # = 1-exp(-theta)
	       x <- Ie/h # = (1-exp(-theta*u)) / (1-exp(-theta))
	       if(log) log(d)+(d-1)*log(x)+log((1-h*x)/(1-h*x^d)) else
	       d*x^(d-1)*(1-h*x)/(1-h*x^d)
	   },
	   "MMH" = { # Martin's version of "MH"
	       h <- -expm1(-theta) # =	1 - exp(-theta)
	       x <- Ie/h #-> h*x = Ie
	       r <- ## log( (1-h*x)/(1-h*x^d) ); 1-h*x = 1-Ie = exp(-u * theta) = exp(-ut)
		   -ut - log1p(-h*x^d)
	       if(log) log(d)+(d-1)*log(x) + r else d*x^(d-1)*exp(r)
	   },
	   stop("impossible method: ", method,". Should never happen"))
}


### Gumbel #####################################################################

### compute the coefficients for polyG

##' Compute the coefficients a_{dk}(theta) involved in the generator derivatives
##' and the copula density of a Gumbel copula
##'
##' @title Coefficients of the polynomial involved in the generator derivatives
##'        and density for Gumbel
##' @param d number of coefficients, d >= 1
##' @param alpha parameter (1/theta) in (0,1];
##'     you maye use  mpfr(alph, precBits = <n_prec>)  for higher precision "Rmpfr*" methods
##' @param method a character string, one of
##'    "sort":          compute coefficients via exp(log()) pulling out the maximum, and sort
##'    "horner":        uses polynEval()
##'    "direct":        brute force approach
##'    "dsSib.<FOO>":   uses dsumSibuya(..., method= "<FOO>")
##' @param log boolean which determines if the logarithm is returned
##' @param verbose logical for method == sort
##' @return a_k(theta, d) = (-1)^{d-k}\sum_{j=k}^d alpha^j * s(d,j) * S(j,k), k in
##'         {1,..,d}
##' @author Marius Hofert and Martin Maechler
##' Note: - This function is known to cause numerical problems, e.g., for d=100,
##'         alpha=0.8
##'       - sign(s(n,k)) = (-1)^{n-k}
coeffG <- function(d, alpha,
                   method = c("sort", "horner", "direct", "dsumSibuya",
                              paste("dsSib", eval(formals(dsumSibuya)$method), sep=".")),
		   log = FALSE, verbose = FALSE)
{
    stopifnot(is.numeric(d), length(d) == 1, d >= 1, length(alpha) == 1,
	      0 < alpha, alpha <= 1)
    a <- numeric(d) # for the a_k(theta, d)'s
    method <- match.arg(method)
    if(method == "dsumSibuya") {
        message("method 'dsumSibuya' is deprecated.\n  ",
        	"use method = 'dsSib.log' instead")
        method <- "dsSib.log"
    }
    Meth <-
	if(grepl("^dsSib", method)) {
	    meth.dsSib <- sub("^dsSib\\.(.*)", "\\1", method)
	    "dsSib"
	} else method
    switch(Meth,
	   "sort" = {
	       ls <- log(abs(Stirling1.all(d))) # log(|s(d, i)|), i=1:d
	       lS <- lapply(1:d, function(n) log(Stirling2.all(n)))
	       ##-> lS[[n]][i] contains log(S(n,i)), i = 1,..,n
               wrong.sign <- integer()
	       for(k in 1:d) { # deal with a_k(theta, d)
		   j <- k:d
		   ## compute b_j = log(alpha^j*|s(d,j)|*S(j,k)), j = k,..,d
		   b <- j * log(alpha) + ls[j] +
		       unlist(lapply(j, function(i) lS[[i]][k]))
		   b.max <- max(b) # max_{j=k,..,d} b_j
		   exponents <- b - b.max # exponents
		   ## compute critical sum (known to be positive)
		   exps <- exp(exponents) # (d-k+1)-many
		   even <- if(k == d) NULL else seq(2, d-k+1, by=2)
		   odd <- seq(1, d-k+1, by=2)
		   sum.neg <- sum(sort(exps[even]))
		   sum.pos <- sum(sort(exps[odd]))
		   sum. <- sum.pos - sum.neg
		   a[k] <- if(log) b.max + log(sum.) else exp(b.max)*sum.
		   if(sum.neg > sum.pos) {
		       if(verbose) message("sum.neg > sum.pos for k = ", k)
		       wrong.sign <- c(wrong.sign, k)
		   }
	       }
	       if(length(wrong.sign))
		   attr(a, "wrong.signs") <- wrong.sign
	       a
	   },
	   "dsSib" = {
	       ## coefficients via dsumSibuya
	       ## a_k(theta,d) = d!/k! * dsumSibuya(d, k, alpha)
	       k <- 1:d
	       ck <- ## c_k := d!/k!
		   if(log) c(0,cumsum(log(d:2)))[d:1]
		   else c(1,cumprod(d:2))[d:1]
	       p <- dsumSibuya(d, k, alpha, method = meth.dsSib, log=log)
               ##   ----------              -------------------
	       if(log) p + ck else p * ck
	   },
	   "horner" = {
	       s.abs <- abs(Stirling1.all(d))
	       k <- 1:d
	       S <- lapply(k, Stirling2.all)## S[[n]][i] contains S(n,i), i = 1,...,n
	       pol <- vapply(k, function(k.) {
		   j <- 0:(d-k.)
		   ## compute coefficients (c_k = |s(d,j+k)|*S(j+k,k))
		   c.j <- s.abs[j+k.] * vapply(j, function(i) S[[i+k.]][k.], 1.)
		   polynEval(c.j, -alpha)
	       }, NA_real_)

	       if(log) k*log(alpha) + log(pol) else alpha^k * pol
	   },
	   "direct" = {
	       s <- Stirling1.all(d) # s(d,1), ..., s(d,d)
	       k <- 1:d
	       S <- lapply(k, Stirling2.all)## S[[n]][i] contains S(n,i), i = 1,...,n
	       vapply(k, function(k.) {
		   j <- k.:d
		   ## extract a column of Stirling2 numbers:
		   S. <- vapply(j, function(i) S[[i]][k.], 1.)
		   sm <- sum(alpha^j * s[j]*S.)
		   if(log) log(abs(sm)) else (-1)^(d-k.)*sm
	       }, NA_real_)
	   },
	   stop(gettextf("unsupported method '%s' in coeffG", method), domain=NA)
	   ) ## switch()
} ## coeffG()

### compute the polynomial for Gumbel

##' Compute the polynomial involved in the generator derivatives and the
##' copula density of a Gumbel copula
##'
##' @title Polynomial involved in the generator derivatives and density for Gumbel
##' @param lx = log(x); where x: evaluation point (vector);
##'        e.g., for copGumbel@dacopula, lx = alpha*log(rowSums(iPsi(u)))
##'        where u = (u_1,..,u_d) is the evaluation point of the density of Joe's copula)
##' @param alpha parameter in (0,1]   alpha := 1/theta = 1 - tau
##' @param d number of summands, >= 1
##' @param method a string, one of
##'   "default"         uses a combination of the other methods
##'   "pois.direct"     uses ppois directly
##'   "pois"            uses ppois with pulling out max
##'   "stirling"        uses the representation via Stirling numbers and once horner
##'   "stirling.horner" uses the representation via Stirling numbers and twice horner
##'    .....            all the method from coeffG(), see there
##' @param verboseUsingRmpfr logical indicating if it should be "messaged" when Rmpfr is used
##'   by the default method.
##' @param log boolean which determines if the logarithm is returned
##' @return \sum_{k=1}^d  a_{dk}(\theta)  x ^ k
##'       = \sum_{k=1}^d  a_{dk} *     exp(lx*k)
##'  where a_{dk}(theta)
##'       = (-1)^{d-k}\sum_{j=k}^d \theta^{-j} s(d,j) S(j,k)
##'       = (d!/k!)\sum_{l=1}^k (-1)^{d-l} \binom{k}{l}\binom{\alpha l}{d}
##' @author Marius Hofert and Martin Maechler
polyG <- function(lx, alpha, d, method= c("default", "default2012", "default2011",
                                "pois", "pois.direct", "stirling", "stirling.horner",
                                coeffG.methods),
                  verboseUsingRmpfr = isTRUE(getOption("copula:verboseUsingRmpfr")),
                  log=FALSE)
{
    stopifnot(length(alpha)==1, 0 < alpha, alpha <= 1,
              d == as.integer(d), d >= 1)
    k <- 1:d
    allMeths <- eval(formals()[["method"]])
    method <- match.arg(method, choices = allMeths)
    Meth <- if(method %in% coeffG.methods) "coeffG" else method
    switch(Meth,
	   "default" =, "default2012" =
       {
	   ## "default2012" compiled by Yongsheng Wang (MSc thesis c/o M.Maechler, April 2012)
	   ## it switches to "Rmpfr" when the accuracy would be less than 5 digits
	   meth2012 <- function(d, alpha, lx) {
	       if (d <= 30) "direct"
	       else if (d <= 50) {
		   if (alpha <= 0.8) "direct" else "dsSib.log"
	       }
	       else if (d <= 70) {
		   if (alpha <= 0.7) "direct" else "dsSib.log"
	       }
	       else if (d <= 90) {
		   if (alpha <= 0.5) "direct"
		   else if (alpha >= 0.8) "dsSib.log"
		   else if (lx <= 4.08) "pois"
		   else if (lx >= 5.4) "direct"
		   else "dsSib.Rmpfr"
	       }
	       else if (d <= 120) {
		   if (alpha < 0.003) "sort"
		   else if (alpha <= 0.4) "direct"
		   else if (alpha >= 0.8) "dsSib.log"
		   else if (lx <= 3.55) "pois"
		   else if (lx >= 5.92) "direct"
		   else "dsSib.Rmpfr"
	       }
	       else if (d <= 170) {
		   if (alpha < 0.01) "sort"
		   else	if (alpha <= 0.3) "direct"
		   else if (alpha >= 0.9) "dsSib.log"
		   else if (lx <= 3.55) "pois"
		   else "dsSib.Rmpfr"
	       }
	       else if (d <= 200) {
		   if (lx <= 2.56) "pois"
		   else if (alpha >= 0.9) "dsSib.log"
		   else "dsSib.Rmpfr"
	       }
	       else "dsSib.Rmpfr"
	   }
	   ## Each lx can -- in principle -- ask for another method ... --> split() by method
	   meth.lx <- vapply(lx, function(lx) meth2012(d, alpha, lx), "")
	   if(verboseUsingRmpfr && (lg <- length(grep("Rmpfr$", meth.lx))))
	       message("Default method chose 'Rmpfr' ", if(lg > 1) paste(lg,"times") else "once")
	   i.m <- split(seq_along(lx), factor(meth.lx))
	   r <- lapply(names(i.m), function(meth)
		       polyG(lx[i.m[[meth]]], alpha = alpha, d = d, method = meth, log = log))
	   lx[unlist(i.m, use.names=FALSE)] <- unlist(r, use.names=FALSE)
	   lx
       },

	   "default2011" = ## first "old" default
       {
	   Recall(lx, alpha=alpha, d=d,
                  method = if(d <= 100) {
                      if(alpha <= 0.54) "stirling"
                      else if(alpha <= 0.77) "pois.direct"
                      else "dsSib.log"
                  } else "pois", # slower but more stable, e.g., for d=150
                  log=log)
       },
           "pois" =
       {
           ## build list of b's
           n <- length(lx)
           x <- exp(lx)                                   # e^lx = x
           lppois <- outer(d-k, x, FUN=ppois, log.p=TRUE) # a (d x n)-matrix; log(ppois(d-k, x))
           llx <- k %*% t(lx)           # also a (d x n)-matrix; k*lx
           labsPoch <- vapply(k, function(j) sum(log(abs(alpha*j-(k-1L)))), NA_real_) # log|(alpha*j)_d|, j=1,..,d
           lfac <- lfactorial(k)        # log(j!), j=1,..,d
           ## build matrix of exponents
           lxabs <- llx + lppois + rep(labsPoch - lfac, n) + rep(x, each = d)
           res <- lssum(lxabs, signFF(alpha, k, d), strict=FALSE)
           if(log) res else exp(res)
       },
           "pois.direct" =
       {
           ## build coefficients
           xfree <- lchoose(alpha*k,d) + lfactorial(d) - lfactorial(k)
           x <- exp(lx)
           lppois <- outer(d-k, x, FUN=ppois, log.p=TRUE) # (length(x),d)-matrix
           klx <- lx %*% t(k)
           exponents <- exp(t(x+klx)+lppois+xfree) # (d,length(x))-matrix
           res <- as.vector(signFF(alpha, k, d) %*% exponents)
           if(log) log(res) else res
       },
           "stirling" =
       {
	   ## implementation of \sum_{k=1}^d a_{dk}(\theta) x^k
	   ## = (-1)^{d-1} * x * \sum_{k=1}^d alpha^k * s(d,k) * \sum_{j=1}^k S(k,j) * (-x)^{j-1}
	   ## = (-1)^{d-1} * x * \sum_{k=1}^d alpha^k * s(d,k) * polynEval(...)
	   ## inner function is evaluated via polynEval
	   x <- exp(lx)
	   s <- Stirling1.all(d) # s(d,1), ..., s(d,d)
	   S <- lapply(k, Stirling2.all) # S[[l]][n] contains S(l,n), n = 1,...,l
	   lst <- lapply(k, function(k.) (-1)^(d-1)*x*alpha^k.*s[k.]*polynEval(S[[k.]],-x))
	   res <- rowSums(matrix(unlist(lst), nrow=length(x)))
	   if(log) log(res) else res
       },
           "stirling.horner" =
       {
	   ## implementation of \sum_{k=1}^d a_{dk}(\theta) x^k
	   ## = (-1)^{d-1} * x * \sum_{k=1}^d alpha^k * s(d,k) * \sum_{j=1}^k S(k,j) * (-x)^{j-1}
	   ## = (-1)^{d-1} * x * \sum_{k=1}^d alpha^k * s(d,k) * polynEval(...)
	   ## polynEval is used twice
	   x <- exp(lx)
	   s <- Stirling1.all(d) # s(d,1), ..., s(d,d)
	   S <- lapply(k, Stirling2.all) # S[[l]][n] contains S(l,n), n = 1,...,l
	   len <- length(x)
           poly <- matrix(unlist(lapply(k, function(k.) polynEval(S[[k.]],-x))), nrow=len) # (len,d)-matrix
           res <- (-1)^(d-1)*alpha*x* vapply(1:len, function(i) polynEval(s*poly[i,], alpha), 1.)
           if(log) log(res) else res
           ## the following code was *not* faster
           ## poly <- t(sapply(k, function(k.) polynEval(S[[k.]],-x))) # (d,len(x))-matrix
           ## coeff <- if(length(x)==1) t(s*poly) else s*poly
           ## res <- (-1)^(d-1)*alpha*x*apply(coeff, 2, polynEval, x=alpha)
       },
           "coeffG" = ## <<< all the different 'coeffG' methods --------------------
       {
           ## note: these methods are all known to show numerical deficiencies
           if(d > 220) stop("d > 220 not yet supported") # would need Stirling2.all(d, log=TRUE)
           ## compute the log of the coefficients:
	   l.a.dk <- coeffG(d, alpha, method=method, log = TRUE)
	   ##	     ~~~~~~			     ~~~~~~~~~~
           ## evaluate the sum
           ## for this, create a matrix B with (k,i)-th entry
           ## B[k,i] = log(a_{dk}(theta)) + k * lx[i],
           ##          where k in {1,..,d}, i in {1,..,n} [n = length(lx)]
           logx <- l.a.dk + k %*% t(lx)
           if(log) {
               ## compute log(colSums(exp(B))) stably (no overflow) with the idea of
               ## pulling out the maxima
               lsum(logx)
           } else colSums(exp(logx))
       },
	   stop(gettextf("unsupported method '%s' in polyG", method), domain=NA)
	   ) # end{switch}
}## {polyG}


### Joe ########################################################################

### sampling a Sibuya(alpha) distribution, R version

##' Sample V from a Sibuya(alpha) distribution with cdf F(n) = 1-1/(n*B(n,1-alpha)),
##' n in IN, with Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the
##' algorithm of Hofert (2011), Proposition 3.2. R version.
##'
##' @title Sampling Sibuya(alpha) distributions
##' @param n  sample size
##' @param alpha parameter in (0,1]
##' @return vector of random variates V
##' @author Marius Hofert, Martin Maechler
rSibuyaR <- function(n, alpha) {
    stopifnot((n <- as.integer(n)) >= 0, length(alpha)==1, 0 < alpha, alpha <= 1)
    V <- numeric(n)
    if(n >= 1) {
        if(alpha == 1) {
            V[] <- 1
        } else {
            u <- runif(n)
            ## FIXME(MM): (for alpha not too close to 1): re-express using 1-u
            l1 <- u <= alpha
            V[l1] <- 1
            i2 <- which(!l1)
            Ginv <- ((1-u[i2])*gamma(1-alpha))^(-1/alpha)
            floorGinv <- floor(Ginv)
            l3 <- (1-1/(floorGinv*beta(floorGinv,1-alpha)) < u[i2])
            V[i2[l3]] <- ceiling(Ginv[l3])
            i4 <- which(!l3)
            V[i2[i4]] <- floorGinv[i4]
        }
    }
    V
}

### state-of-the-art: sampling a Sibuya(alpha) distribution, C version

##' Sample V from a Sibuya(alpha) distribution with cdf F(n) = 1-1/(n*B(n,1-alpha)),
##' n in IN, with Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the
##' algorithm of Hofert (2011), Proposition 3.2. C version.
##'
##' @title Efficiently sampling Sibuya(alpha) distributions
##' @param n sample size (has to be numeric, >= 0)
##' @param alpha parameter in (0,1]
##' @return vector of random variates V
##' @author Martin Maechler
rSibuya <- function(n, alpha) {
    .Call(rSibuya_vec_c, n, alpha)
}

##' Probability mass function of a Sibuya(alpha) distribution
##'
##' @title Probability mass function of a Sibuya(alpha) distribution
##' @param x evaluation point [integer]
##' @param alpha parameter alpha
##' @param log boolean which determines if the logarithm is returned
##' @return p_x = choose(alpha, x) * (-1)^(x-1)
##' @author Marius Hofert and Martin Maechler
dSibuya <- function(x, alpha, log=FALSE)
    if(log) lchoose(alpha, x) else abs(choose(alpha, x))

##' Distribution function of a Sibuya(alpha) distribution
##'
##' @title Distribution function of a Sibuya(alpha) distribution
##' @param x evaluation point [integer]
##' @param alpha parameter alpha
##' @param lower.tail if TRUE, probabilities are P[X <= x], otherwise, P[X > x]
##' @param log.p boolean which determines if the logarithm is returned
##' @return F(x) = 1 - (-1)^x * choose(alpha-1, x)
##' @author Marius Hofert and Martin Maechler
pSibuya <- function(x, alpha, lower.tail=TRUE, log.p=FALSE)
{
    ## F(x) = 1 - 1/(x*Beta(x,1-alpha)) = 1 - (x*beta(x, 1-alpha))^(-1)
    if(log.p) {
        if(lower.tail) # log(1 - 1/(x*beta(x, 1-alpha)))
            log1p(-1/(x*beta(x, 1-alpha)))
        else ## log(1/(x*beta(x, 1-alpha))) = - log(x * beta(..)) =
            -log(x) - lbeta(x, 1-alpha)
    } else { ## no log
        xb <- 1/(x*beta(x, 1-alpha))
        if(lower.tail) 1 - xb else xb
    }
}

### state-of-the art: sampling F01Joe, C version

##' Generate a vector of variates V01 ~ F01 with Laplace-Stieltjes transform
##' ((1-(1-exp(-t))^alpha))^V0. Bridge to R. Used, e.g., to draw several variates
##' from rF01Joe.
##'
##' @title Sampling F01 for Joe's family
##' @param V0 vector of random variates from F0
##' @param parameter alpha = theta0/theta1 in (0,1]
##' @param approx largest number of summands before asymptotics is used
##' @return vector of random variates V01
##' @author Marius Hofert
rF01Joe <- function(V0, alpha, approx) {
    .Call(rF01Joe_vec_c, V0, alpha, approx)
}

### wrapper for inner distribution F for Joe

##' Generate a vector of variates V ~ F with Laplace-Stieltjes transform
##' 1-(1-exp(-t))^alpha.
##'
##' @title Sampling F for Joe
##' @param n number of variates from F
##' @param parameter alpha = theta0/theta1 in (0,1]
##' @return vector of random variates V
##' @author Marius Hofert
rFJoe <- function(n, alpha) rSibuya(n, alpha)

### polynomial evaluation for Joe

##' Inner probability mass function for a nested Joe copula, i.e. a Sibuya sum
##' also used in coeffG() -> polyG() for Gumbel's copula
##'
##' @title Inner probability mass function for a nested Joe copula
##' @param x vector (or number) of natural numbers
##' @param n vector (or number) of natural numbers
##' @param alpha parameter in (0,1]
##' @param method method applied
##'        log:      proper log computation based on lssum
##'        direct:   brute-force evaluation of the sum and its log
##'        Rmpfr, RmpfrM:    multi-precision; the latter *return* multi-prec.
##'        diff:     via forward differences
##'        exp.log:  similar to method = "log", but without *proper/intelligent* log
##' @param log boolean which determines if the logarithm is returned
##' @return p_{x,n} = \sum_{j=1}^n choose(n,j)*choose(alpha*j,x)*(-1)^(x-j)
##'         which is a probability mass function in x on IN with generating function
##'         g(z) = (1-(1-z)^alpha)^n
##' @author Marius Hofert and Martin Maechler
##' note: - p_{x,n} = 0 for x < n;  p_{n,n} = alpha^n
##'       - numerically challenging, e.g., dsumSibuya(100, 96, 0.01) < 0 for all methods
##'  o  Called as  dsumSibuya(d, k, alpha, method = meth.dsSib, log=log)
##'     from coeffG()                      -------------------
dsumSibuya <- function(x, n, alpha,
		       method= c("log", "direct", "diff", "exp.log",
				"Rmpfr", "Rmpfr0", "RmpfrM", "Rmpfr0M"),
                       mpfr.ctrl = list(minPrec = 21, fac = 1.25, verbose=TRUE),
		       log=FALSE)
{
    stopifnot(x == round(x), n == round(n), n >= 0, length(alpha) == 1,
	      0 < alpha, alpha <= 1)
    if((l.x <- length(x)) == 0 || (l.n <- length(n)) == 0)
	return(numeric())
    ## "FIXME": from coeffG(), always have  {x = d, n = 1:d} -- can be faster *not* recycling
    if((len <- l.x) != l.n) { ## do recycle to common length
	len <- max(l.x, l.n)
	if(l.x < len) {
	    x. <- x
	    x <- rep(x, length.out = len)
	}
	else ## if(l.n < len)
	    n <- rep(n, length.out = len)
    }
    if(alpha == 1)
	return(x == n)
    method <- if(missing(method) && is(alpha, "mpfr"))
	"Rmpfr" else match.arg(method)
    switch(method,
	   "log" =
       {
	   ## computes *proper* log based on lssum

	   ## determine the matrix of signs of choose(alpha*j,x)*(-1)^(x-j),
	   ## j in {1,..,m} -- which notably do *not* depend on x !
	   signs <- signFF(alpha, seq_len(max(n)), x)
	   ## for one pair of x and n:
	   f.one <- function(x,n) {
	       if(x < n) return(-Inf)	# = log(0)
	       j <- seq_len(n)
	       lxabs <- lchoose(n, j) + lchoose(alpha*j, x)
	       ## *NON*-strict -- otherwise need try() :
	       lssum(as.matrix(lxabs), signs[j], strict=FALSE)
	   }
	   S <- mapply(f.one, x, n, USE.NAMES = FALSE)
	   if(log) S else exp(S)
       },
	   "direct" =
       {
	   ## brute force evaluation of the sum and its log
	   f.one <- function(x,n) {
	       if(x < n) return(0)
	       j <- seq_len(n)
	       sum(choose(n,j)*choose(alpha*j,x)*(-1)^(x-j))
	   }
	   S <- mapply(f.one, x, n, USE.NAMES = FALSE)
	   if(log) log(S) else S
       },
	   "Rmpfr" =,  "Rmpfr0" =,
	   "RmpfrM" =, "Rmpfr0M" =
       {
	   stopifnot(requireNamespace("Rmpfr"))# need package: classes, methods, e.g. as.numeric(), new()
	   ## as "direct" but using high-precision arithmetic, where
	   ## the precision should be set via alpha = mpfr(*, precBits= .)
	   mayRecall <- !grepl("Rmpfr0", method) ## only if not "Rmpfr0(M)"
	   if(mayRecall) {
	       ## FIXME(?): Hmm, when recalling, alpha becomes more and more precise, even going from n[i] to n[i+1]
	       ## that's not ok.. strictly should only recall for those (x,n) pairs where it's needed
	       stopifnot(is.numeric(minPrec <- mpfr.ctrl$minPrec),
			 is.numeric(fac.pr <- mpfr.ctrl$fac), fac.pr > 1,
			 length(mpfr.ctrl$verbose) == 1)
	   }
	   mayRecall <- !grepl("Rmpfr0", method) ## only if not "Rmpfr0(M)"
	   if(mayRecall) {
	       stopifnot(is.numeric(minPrec <- mpfr.ctrl$minPrec),
			 is.numeric(fac.pr <- mpfr.ctrl$fac), fac.pr > 1,
			 length(mpfr.ctrl$verbose) == 1)
	   }
	   mpfr.0 <- Rmpfr::mpfr(0, precBits = if(!is(alpha, "mpfr")) 64
					       else Rmpfr::getPrec(alpha))
	   ## FIXME: --- speedup possible! ---
	   if(FALSE && l.x == 1 && l.n == x. && all(n == seq_len(x.))) { ## Special case -- from coeffG()
	       message("fast special case ..") ## <- just for now
	       ## change notation (for now)
	       j <- n  # == 1:d
	       n <- x. # == d
	       c.n <- Rmpfr::chooseMpfr.all(n)
	       ca.j <- Rmpfr::chooseMpfr(alpha*j,n)*(-1)^(n-j)
	       f.sp <- function(j) {
		   j. <- seq_len(j)
		   sum(c.n[j.] * ca.j[j.])
	       }
	       stop("fast special case -- is still wrong !")
	       S <- new("mpfr", unlist(lapply(j, f.sp)))
	   } else { ## usual case
	       ## satisfying codetools package checks:
	       mpfr <- Rmpfr::mpfr; roundMpfr <- Rmpfr::roundMpfr
	       chooseMpfr <- Rmpfr::chooseMpfr; chooseMpfr.all <- Rmpfr::chooseMpfr.all
	       f.1 <- function(x,n, alpha) {
		   if(x < n) return(mpfr.0)
                   pr. <- .dsSib.mpfr.prec(x, n, alpha)
                   ##     ================
                   alpha <- mpfr(alpha, precBits = pr.)
		   j <- seq_len(n)
		   ##if(TRUE) { # "old"/for debug:
		   ## now do this in parts {to analyze}:
		   ## sum(chooseMpfr.all(n)*chooseMpfr(alpha*j,x)*(-1)^(x-j))
		   c.n <- chooseMpfr.all(n)
		   ca.j <- chooseMpfr(alpha*j,x) * (-1)^(x-j)
		   trms <- c.n * ca.j
		   S <- sum(trms) # <-- sum of *huge* (partly alternating) terms getting almost 0
		   ## } else {	      # "new", typically faster
### FIXME: Use *faster*
		   ##	  ## sumBinomMpfr(n, FUN, n0, alternating=TRUE, precBits) is not yet available
		   ##	  S <- sum(chooseZ(n,j) * (-1)^(x-j) * chooseMpfr(alpha * j, x))
		   ## }
		   if(mayRecall) {
		       recall <- (S <= 0)
		       if(recall) { ## complete loss of precision -- recall with higher precision
			   MSG <- " |--> S < 0"
			   newPrec <- round(fac.pr * pr.)
		       } else {		# S > 0 :
			   bitLoss <- round(as.numeric(log2(max(abs(trms)) / S )), 1)
			   recall <- (pr. - bitLoss < minPrec) ## less than 'minPrec' bits precision left
			   if(recall)
			       MSG <- paste("-> bit loss =", bitLoss)
			   newPrec <- round(max(minPrec + bitLoss, fac.pr * pr.))
			   ## newPrec <- round(fac.pr * pr.)
		       }
		       if(recall) {
			   if(mpfr.ctrl$verbose)
                               message(sprintf("dsumSibuya(%d, %d, alpha= %g [pr = %d], method='%s'): %s ==> recalled w/ new prec. %d",
					   x, n, as.numeric(alpha), as.integer(pr.), method, MSG, newPrec))
			   dsumSibuya(x,n, alpha = roundMpfr(alpha, newPrec), mpfr.ctrl=mpfr.ctrl,
				      method="RmpfrM", log=FALSE)
		       } else S
		   } else S
	       }## end{ f.1() }
	       S <- new("mpfr", mapply(f.1, x, n, alpha, USE.NAMES=FALSE))
	   }
	   if(grepl("M$", method)) ## "RmpfrM" or "Rmpfr0M"
	       if(log) log(S) else S
	   else
	       as.numeric(if(log) log(S) else S)
       },
	   "diff" =
       {
	   S. <- mapply(function(x,n) diff(choose(n:0*alpha, x), differences=n) * (-1)^x,
			x, n, USE.NAMES = FALSE)
	   if(log) log(S.) else S.
       },
	   "exp.log" =
       {
	   ## similar to method = "log", but without *proper/intelligent* log
	   ## and inefficient due to the signs (old version)
	   f.one <- function(x,n) {
	       if(x < n) return(0)
	       j <- seq_len(n) ## indices of the summands
	       signs <- (-1)^(j+x)
	       ## determine the signs of choose(j*alpha,x) for each component of j
	       to.subtract <- 0:(x-1)
	       sig.choose <-
		   vapply(j, function(l) prod(sign(l*alpha-to.subtract)),
			  1., USE.NAMES=FALSE)
	       signs <- signs*sig.choose
	       binom.coeffs <- exp(lchoose(n,j) + lchoose(j*alpha,x))
	       sum(signs*binom.coeffs)
	   }
	   S <- mapply(f.one, x, n, USE.NAMES = FALSE)
	   if(log) log(S) else S
       },
	   ## otherwise
	   stop(gettextf("unsupported method '%s' in dsumSibuya", method)), domain=NA)
}

..dsSib.mpfr.prec <- function(d,k,alpha)
{
    ## no checking here on purpose -- note that this works *vectorized*
    ldk <- log(d - k + 1)
    la <- log(alpha)
    sa <- sqrt(alpha)
    ## .dsSib.mpfr.precEXPR  from above:
    yhat <- 0.606778486870861  + 1.00415795049476 * log(k) +
        -0.0115467115092516 * ldk + -0.236630339094924 * la + 6.84638655885601 * sa +
            -7.59637383430576 * alpha + -0.529181206860549 * ldk * sa +
                0.579077302168194 * ldk * alpha + 4.07566657020875 * la * alpha
    pmax(64, exp(yhat + 0.18))
}

.dsSib.mpfr.prec <- function(d, k, alpha)
{
    stopifnot(is.numeric(d), is.numeric(k), is.numeric(alpha),
              length(d) == 1, length(k) == 1, length(alpha) == 1,
              d == as.integer(d), k == as.integer(k),
              d >= k, 0 < alpha, alpha < 1)
    if(alpha < 0.001)
	stop("mpfr precision needed for	 alpha < 0.001	is not yet available.\n",
	     " Please report to maintainer(\"copula\") if you need this to be changed")
    ## @MM: --> ~/R/MM/Pkg-ex/copula/dsSibMpfr-prec/dsSibMpfr-ex.R
    p <-
        if (k <= 5)
            64
        else if (alpha >= 0.2) {
            if (d-k >= 170)
                64
            else if (d-k >= 90) {
                if (d <= 120)
                    64
            } else if (d-k > 2) {
                if ((d - 70*alpha) < 9)
                    64
            }
        }
    if(is.null(p)) ..dsSib.mpfr.prec(d,k,alpha) else p
}


### polynomial evaluation for Joe

##' Compute the polynomial involved in the generator derivatives and the
##' copula density of a Joe copula
##'
##' @title Polynomial involved in the generator derivatives and density for Joe
##' @param lx (log) evaluation point (lx is meant to be log(x) for some x which
##'        was used earlier; e.g., for copJoe@dacopula, lx = log(h(u)/(1-h(u))) for
##'        h(u) = \prod_{j=1}^d(1-(1-u_j)^theta), where u = (u_1,..,u_d) is the
##'        evaluation point of the density of Joe's copula)
##' @param alpha parameter alpha ( := 1/theta ) in (0,1]
##' @param d number of summands
##' @param method different methods, can be
##'        "log.poly" intelligent log version
##'        "log1p"    additonally uses log1p
##'        "poly"     brute force log version
##' @param log boolean which determines if the logarithm is returned
##' @return \sum_{k=1}^d a_{d,k}(theta) exp((k-1)*lx) = \sum_{k=0}^{d-1} a_{d,k+1}(theta) x^k
##'         where a_{d,k}(theta) = S(d,k)*(k-1-alpha)_{k-1} = S(d,k)*Gamma((1:d)-alpha)/Gamma(1-alpha)
##'         Note: these a_{d,k}(theta) here are not those of Hofert, Maechler, McNeil (2012)
##'               (they are the a_{d,k-1}(theta))
##' @author Marius Hofert and Martin Maechler
polyJ <- function(lx, alpha, d, method=c("log.poly","log1p","poly"), log=FALSE) {
    stopifnot(length(alpha)==1, 0 < alpha, alpha <= 1)
    if(!length(lx)) return(numeric())
    ## compute the log of the coefficients a_{dk}(theta)
    if(d > 220) stop("d > 220 not yet supported")# would need Stirling2.all(d, log=TRUE)
    k <- 1:d
    l.a.k <- log(Stirling2.all(d)) + lgamma(k-alpha) - lgamma(1-alpha) # log(a_{dk}(theta)), k = 1,..,d
    ## evaluate polynomial via exp( log(<poly>) )
    ## for this, create a matrix B with (k,i)-th entry B[k,i] = log(a_{dk}(theta)) + (k-1) * lx[i],
    ## where k in {1,..,d}, i in {1,..,n} [n = length(lx)]
    B <- l.a.k + (k-1) %*% t(lx)
    method <- match.arg(method)
    switch(method,
           "log.poly" = {
               ## stably compute log(colSums(exp(B))) (no overflow)
               ## Idea:
               ## (1) let b_k := log(a_{dk}(theta)) + (k-1)*lx and b_{max} := argmax{b_k}.
               ## (2) \sum_{k=1}^d a_{dk}(theta)\exp((k-1)*lx) = \sum_{k=1}^d \exp(log(a_{dk}(theta))
               ##     + (k-1)*lx) = \sum_{k=1}^d \exp(b_k) = \exp(b_{max})*\sum_{k=1}^d
               ##     \exp(b_k-b_{max})
               ## (3) => log(\sum...) = b_{max} + log(\sum_{k=1}^d \exp(b_k-b_{max}))
               if(log) lsum(B) else exp(lsum(B))
           },
           "log1p" = {
               ## use log(1 + sum(<smaller>)) = log1p(sum(<smaller>)),
               ## but we don't expect it to make a difference
               im <- apply(B, 2, which.max) # indices (vector) of maxima
               n <- length(lx) ; d1 <- d-1L
               max.B <- B[cbind(im, seq_len(n))] # get max(B[,i])_{i=1,..,n} == apply(B, 2, max)
               B.wo.max <- matrix(B[unlist(lapply(im, function(j) k[-j])) +
                                    d*rep(0:(n-1), each = d1)], d1, n) # matrix B without maxima
               res <- max.B + log1p(colSums(exp(B.wo.max - rep(max.B, each = d1))))
               if(log) res else exp(res)
           },
           "poly" = {
               ## brute force ansatz
               res <- colSums(exp(B))
               if(log) log(res) else res
           },
       stop(gettextf("unsupported method '%s' in polyJ", method)), domain=NA)
}

##' Circular/Rational function	(1 - x^d)/(1 - x) for x ~~ 1, i.e.,
##' compute (1 - x^d)/(1 - x) = (1 - (1-e)^d) / e   for	 e = 1-x (<< 1) and integer d
##'
##' @title Circular/Rational function  (1 - (1-e)^d) / e  {incl. limit e -> 0}
##' @param e numeric vector in [0, 1]
##' @param d integer (scalar), >= 1
##' @return (1 - (1-e)^d) / e
##' @author Martin Maechler, Date: 25 Jul 2011
circRat <- function(e, d)
{
### TODO (?):  improve "log=TRUE", for e ~= 1:	log(circRat(e, d)) = log1p(-x^d) - log(e)
    stopifnot(length(d) == 1, d == as.integer(d), d >= 1)
    if(d <= 6)
	switch(d,
	       1-0*e, ## <<- d = 1
	       2-e,   ## <<- d = 2: 1 2 1
	       3-e*(3-e),#   d = 3: 1 3 3 1
	       4-e*(6-e*(4-e)),		       ## d = 4: 1 4 6 4 1
	       5-e*(10-e*(10-e*(5-e))),	       ## d = 5: 1 5 10 10 5 1
	       6-e*(15-e*(20-e*(15-e*(6-e))))  ## d = 6: 1 6 15 20 15 6 1
	       )
    else { ## d >= 7 ---------------
	r <- e
	eps <- .Machine$double.eps
	if(any(l1 <- ((d1 <- (d-1)/2)*e < eps)))
	    r[l1] <- d
	if(any(l2 <- !l1 & ((d2 <- (d-2)/3)*(e2 <- e*e) < eps)))
	    r[l2] <- d*(1 - d1*e[l2])
	if(any(l3 <- !l1 & !l2 & ((d3 <- (d-3)/4)*e*e2 < eps)))
	    r[l3] <- d*(1 - d1*e[l3]*(1 - d2*e[l3]))
	## and for the remaining ones, we afford a little precision loss:
	if(any(lrg <- !l1 & !l2 & !l3)) {
	    e <- e[lrg]
	    r[lrg] <- (1 - (1-e)^d)/e
	}
	r
    }
}

##' @title tau(theta) for Joe's copula
##' @param theta (vector of) copula parameters in [-1,1] (in [0,1] for tau >= 0)
##' @param method string specifying the computational method.
##' @param noTerms number of terms to use for method = "sum".
##' @return vector of same length as theta with   tau(theta[.])
##' @author Martin Maechler
tauJoe <- function(theta, method = c("hybrid", "digamma", "sum"), noTerms=446)
{
    method <- match.arg(method)
    switch(method,
	   "hybrid" = {
	       digam1 <- digamma(1) ## == - (Euler's) gamma = -0.5772157
	       trigam1 <- pi^2/6 ##  == psigamma(1, d=1) = trigamma(1)
	       vapply(theta, function(th) {
		   if(th == 2) return(2 - trigam1)
		   if(th > 1e17) return(1)
		   q <- 2/th
		   tol <- 1.5e-5 ## MM: from experimentation (Lnx, 64-bit)
				 ## ---> ~/R/MM/Pkg-ex/copula/tauJoe.R
		   dt <- if(abs(e <- q-1) < tol)## th ~= 2: |2-th| < tol*th
		       -(digam1 + q*trigam1 + e/2*psigamma(1, deriv=2))
		   else
		       (digam1 - q*digamma(q))/e
		   1 + q*(dt + digamma(1+q))
	       }, 0.)
	   },
	   "digamma" = {
	       digam1 <- digamma(1) ## == - (Euler's) gamma = -0.5772157
	       vapply(theta, function(th) {
		   if(th == 2) return(2 - pi^2/6)
		   q <- 2/th
		   ## A <- 1/(q*(q-1)) ## only for q != 1,	 i.e.,	th != 2
		   ## 1 + 4/th^2 * (A*digam1 + digamma(1+q)/q + digamma(q)/(1-q))
		   ##  q^2 = 4/th^2; q*A = 1/(q-1)
		   ## 1 + q * (digam1/(q-1) + digamma(1+q) + q/(1-q)*digamma(q))
		   1 + q*((digam1 - q*digamma(q))/(q-1) + digamma(1+q))
	       }, 0.)
	   },
	   "sum" = {
	       k <- noTerms:1
	       sapply(theta,
		      function(th) {
			  tk2 <- th*k + 2
			  1 - 4*sum(1/(k*tk2*(tk2 - th)))
			  ## ==... (1/(k*(th*k+2)*(th*(k-1)+2)))
		      })
	   },
	   stop("unsupported method: ", method))
}

### Misc #######################################################################

## Determine the \dQuote{implied} copula dimension from \code{u}.
## in the same sense that many copula functions use
##  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)


##' @title Implied copula dim()ension
##' @param u
##' @return integer
##' @author Martin Maechler
dimU <- function(u) {
    if(!is.null(d <- dim(u))) d[2L] else length(u)
}



##' Conditional copula function C(u[,d]|u[,1],...,u[,d-1])
##'
##' @title Conditional copula function
##' @param u (n x d)-matrix of evaluation points (first d-1 columns are conditioned on)
##' @param cop an outer_nacopula
##' @param n.MC Monte Carlo sample size
##' @param log if TRUE the logarithm of the conditional copula is returned
##' @author Marius Hofert
cacopula <- function(u, cop, n.MC=0, log=FALSE) {
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are supported")

    .Deprecated("cCopula")
    d <- ncol(u)
    drop(rtrafo(u, cop=cop, j.ind = d, n.MC=n.MC, log=log))
}

##' Function which computes absdPsi via Monte Carlo
##'
##' @title Computing the absolute value of the generator derivatives via Monte Carlo
##' @param t evaluation points
##' @param family Archimedean family (name or object)
##' @param theta parameter value
##' @param degree order of derivative
##' @param n.MC Monte Carlo sample size
##' @param method different methods
##'        log:         proper log using lsum
##'        direct:      direct evaluation of the sum
##'        pois.direct: directly uses the Poisson density
##'        pois:        intelligently uses the Poisson density with lsum
##' @param log if TRUE the log of absdPsi is returned
##' @param is.log.t if TRUE, the argument t contains log(<mathematical t>)
##' @author Marius Hofert, Martin Maechler
##' Note: absdPsiMC(0) is always finite, although, theoretically, absdPsi(0) may
##'       be Inf (e.g., for Gumbel and Joe)
absdPsiMC <- function(t, family, theta, degree=1, n.MC,
                      method=c("log", "direct", "pois.direct", "pois"),
                      log = FALSE, is.log.t = FALSE)
{
    res <- numeric(length(t))
    V <- getAcop(family)@V0(n.MC, theta)
    Vt <- if(is.log.t) function(tt) V %*% t(exp(tt)) else function(tt) V %*% t(tt)
    method <- match.arg(method)
    switch(method,
	   ## the following is not faster than "log":
	   ## "default" = { # basically, use "direct" if numerically not critical and "log" otherwise
	   ##                lx <- -V %*% t(t) + degree*log(V)
	   ##                explx <- exp(lx) # (n.MC, n)-matrix containing the summands
	   ##                explx0 <- explx==0 # can be TRUE due to t == Inf or t finite but too large
	   ##                t.too.large <- unlist(lapply(1:n, function(x) any(explx0))) # boolean vector of length n indicating which column of explx contains zeros
	   ##                r1 <- colMeans(explx[,!t.too.large, drop=FALSE])
	   ##                res[!t.too.large] <- if(log) log(r1) else r1
	   ##                r2 <- lsum(lx[,t.too.large, drop=FALSE] - log(n.MC))
	   ##                res[t.too.large] <- if(log) r2 else exp(r2)
	   ##                res[is.infinite(t)] <- if(log) -Inf else 0
	   ##                res
	   ##            },
           "log" = { # intelligent log
               iInf <- is.infinite(t)
               res[iInf] <- -Inf # log(0)
               if(any(!iInf))
                   res[!iInf] <- lsum(-Vt(t[!iInf]) + degree*log(V) - log(n.MC))
               if(log) res else exp(res)
           },
           "direct" = { # direct method
               lx <- -Vt(t) + degree*log(V)
               res <- colMeans(exp(lx)) # can be all zeros if lx is too small [e.g., if t is too large]
               if(log) log(res) else res
           },
           "pois.direct" = {
               m.poi <- colMeans(dpois(degree, lambda=Vt(t)))
               ## is.log.t:  "t^degree" = exp(t)^degree = exp(t * degree)
               res <- factorial(degree)*(if(is.log.t) exp(-t * degree) else t^-degree) * m.poi
               if(log) log(res) else res
           },
           "pois" = {
               iInf <- is.infinite(t)
               res[iInf] <- -Inf # log(0)
               if(any(!iInf)) {
                   t <- t[!iInf]
                   lpoi <- dpois(degree, lambda=Vt(t), log=TRUE) # (n.MC, length(t))-matrix
                   b <- -log(n.MC) + lfactorial(degree) - degree*rep(if(is.log.t)t else log(t), each=n.MC) + lpoi
                   res[!iInf] <- lsum(b)
               }
               if(log) res else exp(res)
           },
	   stop(gettextf("unsupported method '%s' in absdPsiMC", method)), domain=NA)
}

psiDabsMC <- function(t, family, theta, degree=1, n.MC,
                      method=c("log", "direct", "pois.direct", "pois"),
                      log = FALSE, is.log.t = FALSE)
{
    .Deprecated("absdPsiMC")
    absdPsiMC(t, family=family, theta=theta, degree=degree, n.MC=n.MC,
                      method=method, log=log, is.log.t)
}

### Non-numerics ###############################################################


setGeneric("setTheta", function(x, value, ...) standardGeneric("setTheta"))

##' @title Setting the parameter in a copula
##' @param x acopula
##' @param value parameter value
##' @param na.ok logical indicating if NA values are ok for theta
##' @param noCheck logical indicating if parameter constraints should be checked
##' @return acopula with theta set to value
##' @author Martin Maechler
setMethod("setTheta", "acopula",
	  function(x, value, na.ok = TRUE, noCheck = FALSE, ...)
      {
	  stopifnot(is.numeric(value) | (ina <- is.na(value)))
	  if(ina) {
	      if(!na.ok) stop("NA value, but 'na.ok' is not TRUE")
	      value <- NA_real_
	  }
	  if(ina || noCheck || x@paraConstr(value)) ## parameter constraints are fulfilled
	      x@theta <- value
	  else
	      stop("theta (=", format(value), ") does not fulfill paraConstr()")
	  x
      })
setMethod("setTheta", signature(x="outer_nacopula", value="numeric"),
	  function(x, value, na.ok = TRUE, noCheck = FALSE, ...) {
	      x@copula <- setTheta(x@copula, value, na.ok=na.ok, noCheck=noCheck)
	      x
	  })

## TODO: setTheta - using a list of thetas
## setMethod("setTheta", signature(x="outer_nacopula", value="list"),
## 	  function(x, value, na.ok = TRUE, noCheck = FALSE) {
##               ... x@copula <- setTheta(x@copula, value, na.ok=na.ok, noCheck=noCheck)
##           })

setMethod("setTheta", "copula",
	  function(x, value, na.ok = TRUE, noCheck = FALSE, ...)
      {
	  stopifnot(is.numeric(value) | (ina <- is.na(value)))
	  if(any(ina)) {
	      if(!na.ok) stop("NA value, but 'na.ok' is not TRUE")
	      ## vectorized (and partial)  value <- NA_real_
	      if(!is.double(value)) storage.mode(value) <- "double"
	  }
	  ## not using has.par.df(), as have "ellipCopula" method below
	  if(is(x, "tevCopula")) {
	      p <- (!x@df.fixed) + 1
	      if(length(value) != p)
		  stop(gettextf("'length(value)' must be %d for tevCopula(dim=%d)",
				p, x@dimension), domain=NA)
	  }
	  ##if(ina || noCheck || x@paraConstr(value)) ## parameter constraints are fulfilled
	  if(all(ina) || noCheck || {
	      all(is.na(value) | (x@param.lowbnd <= value & value <= x@param.upbnd))
	  }) ## parameter constraints are fulfilled
	      x@parameters[seq_along(value)] <- value
	  else
	      stop(gettextf("theta (=%s) is not inside parameter bounds",
			    format(value)), domain=NA)
	  x
      })


## NB:  for tCopula(df.fixed = FALSE), value now is (rho,df)
setMethod("setTheta", "ellipCopula",
	  function(x, value, na.ok = TRUE, noCheck = FALSE, ...)
      {
	  stopifnot(is.numeric(value) | (ina <- is.na(value)))
	  if(any(ina)) {
	      if(!na.ok) stop("NA value, but 'na.ok' is not TRUE")
	      ## vectorized (and partial)  value <- NA_real_
	      if(!is.double(value)) storage.mode(value) <- "double"
	  }
          df.f <- if(is(x, "tCopula")) x@df.fixed else TRUE
          p <- npar.ellip(x@dimension, dispstr = x@dispstr, df.fixed = df.f)
          if(length(value) != p)
              stop(gettextf("'length(value)' must be %d for this elliptical copula (dim=%d, dispstr=\"%s\")",
                            p, x@dimension, x@dispstr), domain=NA)
	  if(all(ina) || noCheck || {
	      all(is.na(value) | (x@param.lowbnd <= value & value <= x@param.upbnd))
	  }) ## parameter constraints are fulfilled
	      x@parameters[seq_along(value)] <- value
	  else
	      stop(gettextf("theta (=%s) is not inside parameter bounds",
			    format(value)), domain=NA)
	  x
      })


##' Construct "paraConstr" function from an "interval"
##'
##' @title Construct "paraConstr" function from an "interval"
##' @param int interval
##' @return parameter constraint function
##' @author Martin Maechler
mkParaConstr <- function(int) {
    stopifnot(is(int, "interval")) # for now
    is.o <- int@open
    eL <- substitute(LL <= theta, list(LL = int[1])); if(is.o[1]) eL[[1]] <-
        as.symbol("<")
    eR <- substitute(theta <= RR, list(RR = int[2])); if(is.o[2]) eR[[1]] <-
        as.symbol("<")
    bod <- substitute(length(theta) == 1 && LEFT && RIGHT,
                      list(LEFT = eL, RIGHT= eR))
    as.function(c(alist(theta=), bod), parent.env(environment()))
    ## which is a fast version of
    ## r <- function(theta) {}
    ## environment(r) <- parent.env(environment())
    ## body(r) <- bod
    ## r
}

printAcopula <- function(x, slots = TRUE, indent = 0,
                         digits = getOption("digits"), width = getOption("width"), ...)
{
    cl <- class(x)
    cld <- getClassDef(cl)
    stopifnot(indent >= 0, extends(cld, "acopula"))
    ch.thet <- {
        if(!all(is.na(x@theta)))## show theta
	    paste0(", theta= (",
		   paste(sapply(x@theta, format, digits=digits), collapse=", "), ")")
        else ""
    }
    bl <- paste(rep.int(" ",indent), collapse="")
    cat(sprintf('%sArchimedean copula ("%s"), family "%s"%s\n',
                bl, cl, x@name, ch.thet))
    if(slots) {
        nms <- slotNames(cld)
        nms <- nms[!(nms %in% c("name", "theta"))]
        i2 <- indent+2
        cat(bl, " It contains further slots, named\n",
            paste(strwrap(paste(dQuote(nms),collapse=", "),
                          width = 0.95 * (width-2), indent=i2, exdent=i2),
                  collapse="\n"), "\n",
            sep="")
    }
    invisible(x)
}
setMethod(show, "acopula", function(object) printAcopula(object))

## This is now exported & has help file --> ../man/printNacopula.Rd :
printNacopula <-
    function(x, labelKids = NA, deltaInd = if(identical(labelKids,FALSE)) 5 else 3,
             indent.str="",
             digits = getOption("digits"), width = getOption("width"), ...)
{
    cl <- class(x)
    stopifnot(deltaInd >= 0, is.character(indent.str), length(indent.str) == 1,
              extends(cl, "nacopula"))
    mkBlanks <- function(n) paste(rep.int(" ", n), collapse="")
    bl <- mkBlanks(nIS <- nchar(indent.str))

    ccl <- if(extends(cl, "outer_nacopula"))
        paste0('"',cl,'" of dim. ', dim(x)) else paste0('"',cl,'"')
    ## cat(sprintf(" __deltaInd = %d, nIS = %d__ ", deltaInd, nIS))
    ch1 <- sprintf("%sNested Archimedean copula (%s), with ",
                   indent.str, ccl)
    ch2 <- if(length(c.j <- x@comp)) {
        sprintf("slot \n%s'comp'   = %s", bl,
		paste0("(",paste(c.j, collapse=", "),")"))
    } else "empty slot 'comp'"
    cat(ch1, ch2, sprintf("  and root\n%s'copula' = ", bl), sep="")
    printAcopula(x@copula, slots=FALSE, digits=digits, width=width, ...)
    nk <- length(kids <- x@childCops)
    if(nk) {
        cat(sprintf("%sand %d child copula%s\n", bl, nk, if(nk > 1)"s" else ""))
        doLab <- if(is.na(labelKids)) nk > 1 else as.logical(labelKids)
        if(doLab) {
            hasNms <- !is.null(nms <- names(kids))
            lab <- if(hasNms) paste0(nms,": ") else paste0(seq_len(nk),") ")
        }
        bl <- mkBlanks(nIS + deltaInd)
        for(ic in seq_along(kids))
            printNacopula(kids[[ic]], deltaInd=deltaInd,
                          indent.str = paste0(bl, if(doLab) lab[ic]),
                          labelKids=labelKids, digits=digits, width=width, ...)
    }
    else
        cat(sprintf("%sand *no* child copulas\n", bl))
    invisible(x)
}

setMethod(show, "nacopula", function(object) printNacopula(object))

##' Get one of our "acopula" family objects by name
##'
##' @title Get one of our "acopula" family objects by name
##' @param family either character string (short or longer form of
##'	 copula family name), an "acopula" family object,
##'      *or* an object inheriting from "archmCopula"
##' @param check logical indicating if the class of the return value should
##' be checked.
##' @return one of our "acopula" objects
##' @author Martin Maechler
getAcop <- function(family, check=TRUE) {
    if(is.character(family)) {
	stopifnot(length(family) == 1)
	if((nf <- nchar(family)) <= 2) # it's a short name
	    family <- .ac.longNames[family]
	else if (nf >= 8 && grepl("Copula$", family))
	    family <- names(which(.ac.classNames == family))
	COP <- get(.ac.objNames[family]) # envir = "package:copula"
	if(check && !is(COP, "acopula"))
	    stop(paste("invalid acopula-family object, family=",family))
	COP
    } else {
	cl <- getClass(class(family))# so the extends(.) below are fast
	if(extends(cl, "acopula"))
	    family
	else if(extends(cl, "archmCopula")) {
	    if(extends(cl, "indepCopula"))## FIXME? do not want full family object!
		stop("independence copula not implemented as Archimedean family")
	    ## now use short family names
	    getAcop(if(extends(cl, "claytonCopula")) "C" else
		    if(extends(cl, "frankCopula"))   "F" else
		    if(extends(cl, "amhCopula"))     "A" else
		    if(extends(cl, "joeCopula"))     "J" else
		    if(extends(cl, "gumbelCopula"))  "G" else
		    stop("invalid archmCopula class: ", cl))
	}
	else stop("'family' must be an \"archmCopula\" or \"acopula\" object or family name")
    }
}


coeffG.methods <- eval(formals(coeffG)$method)# - namespace hidden
## --> accesses formals(dsumSibuya) .. hence at end
