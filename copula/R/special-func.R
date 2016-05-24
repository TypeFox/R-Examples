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


### Special Numeric Functions -- related to Archimedean Copulas
### -------------------------
### ( were "other numeric utilities" in ./aux-acopula.R )

polynEval <- function(coef, x) .Call(polyn_eval, coef, x)


##' @title Compute  f(a) = log(1 - exp(-a))  stably
##' @param a numeric vector of positive values
##' @param cutoff  log(2) is optimal, see  Maechler (201x) .....
##' @return f(a) == log(1 - exp(-a)) == log1p(-exp(-a)) == log(-expm1(-a))
##' @author Martin Maechler, May 2002 .. Aug. 2011
##' @references Maechler(2012)
##' Accurately Computing log(1 - exp(-|a|)) Assessed by the Rmpfr package.
##' http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
## MM: ~/R/Pkgs/Rmpfr/inst/doc/log1mexp-note.Rnw
##--> ../man/log1mexp.Rd
log1mexp <- function(a, cutoff = log(2)) ## << log(2) is optimal >>
{
    if(has.na <- any(ina <- is.na(a))) {
	y <- a
	a <- a[ok <- !ina]
    }
    if(any(a < 0))## a == 0  -->  -Inf	(in both cases)
	warning("'a' >= 0 needed")
    tst <- a <= cutoff
    r <- a
    r[ tst] <- log(-expm1(-a[ tst]))
    r[!tst] <- log1p(-exp(-a[!tst]))
    if(has.na) { y[ok] <- r ; y } else r
}

##' @title Compute  f(x) = log(1 + exp(x))  stably and quickly
##--> ../man/log1mexp.Rd
log1pexp <- function(x, c0 = -37, c1 = 18, c2 = 33.3)
{
    if(has.na <- any(ina <- is.na(x))) {
	y <- x
	x <- x[ok <- !ina]
    }
    r <- exp(x)
    if(any(i <- c0 < x & (i1 <- x <= c1)))
	r[i] <- log1p(r[i])
    if(any(i <- !i1 & (i2 <- x <= c2)))
	r[i] <- x[i] + 1/r[i] # 1/exp(x) = exp(-x)
    if(any(i3 <- !i2))
	r[i3] <- x[i3]
    if(has.na) { y[ok] <- r ; y } else r
}


##' The sign of choose(alpha*j,d)*(-1)^(d-j) vectorized in j
##'
##' @title The sign of choose(alpha*j,d)*(-1)^(d-j)
##' @param alpha alpha (scalar) in (0,1]
##' @param j integer vector in {0,..,d}
##' @param d integer (scalar) >= 0
##' @return sign(choose(alpha*j,d)*(-1)^(d-j))
##' @author Marius Hofert
##' Note: If alpha=1, then
##'       sign( choose(alpha*j, d)*(-1)^(d-j) ) == (-1)^(d-j) if j > d and
##'       sign( choose(alpha*j, d)*(-1)^(d-j) ) == 0 if j = 0
signFF <- function(alpha, j, d) {
    stopifnot(0 < alpha, alpha <= 1, d >= 0, 0 <= j)
    res <- numeric(length(j))
    if(alpha == 1) {
	res[j == d] <- 1
	res[j > d] <- (-1)^(d-j)
    } else {
	res[j > d] <- NA # the formula below does not hold {TODO: find correct sign}
        ## we do not need them in dsumSibuya() and other places...
	x <- alpha*j
	ind <- x != floor(x)
	res[ind] <- (-1)^(j[ind]-ceiling(x[ind]))
    }
    res
}

##' Properly compute log(x_1 + .. + x_n) for a given (n x d)-matrix of n row
##' vectors log(x_1),..,log(x_n) (each of dimension d)
##' Here, x_i > 0  for all i
##' @title Properly compute the logarithm of a sum
##' @param lx (n,d)-matrix containing the row vectors log(x_1),..,log(x_n)
##'        each of dimension d
##' @param l.off the offset to subtract and re-add; ideally in the order of
##'        the maximum of each column
##' @return log(x_1 + .. + x_n) [i.e., OF DIMENSION d!!!] computed via
##'         log(sum(x)) = log(sum(exp(log(x))))
##'         = log(exp(log(x_max))*sum(exp(log(x)-log(x_max))))
##'         = log(x_max) + log(sum(exp(log(x)-log(x_max)))))
##'         = lx.max + log(sum(exp(lx-lx.max)))
##'         => VECTOR OF DIMENSION d
##' @author Marius Hofert, Martin Maechler
lsum <- function(lx, l.off) {
    rx <- length(d <- dim(lx))
    if(mis.off <- missing(l.off)) l.off <- {
	if(rx <= 1L)
	    max(lx)
	else if(rx == 2L)
	    apply(lx, 2L, max)
    }
    if(rx <= 1L) { ## vector
	if(is.finite(l.off))
	    l.off + log(sum(exp(lx - l.off)))
	else if(mis.off || is.na(l.off) || l.off == max(lx))
	    l.off # NA || NaN or all lx == -Inf, or max(.) == Inf
	else
	    stop("'l.off  is infinite but not == max(.)")
    } else if(rx == 2L) { ## matrix
	if(any(x.off <- !is.finite(l.off))) {
	    if(mis.off || isTRUE(all.equal(l.off, apply(lx, 2L, max)))) {
		## we know l.off = colMax(.)
		if(all(x.off)) return(l.off)
		r <- l.off
		iok <- which(!x.off)
		l.of <- l.off[iok]
		r[iok] <- l.of + log(colSums(exp(lx[,iok,drop=FALSE] -
						     rep(l.of, each=d[1]))))
		r
	    } else ## explicitly specified l.off differing from colMax(.)
		stop("'l.off' has non-finite values but differs from default max(.)")
	}
	else
	    l.off + log(colSums(exp(lx - rep(l.off, each=d[1]))))
    } else stop("not yet implemented for arrays of rank >= 3")
}

##' Properly compute log(x_1 + .. + x_n) for a given matrix of column vectors
##' log(|x_1|),.., log(|x_n|) and corresponding signs sign(x_1),.., sign(x_n)
##' Here, x_i is of arbitrary sign
##' @title compute logarithm of a sum with signed large coefficients
##' @param lxabs (d,n)-matrix containing the column vectors log(|x_1|),..,log(|x_n|)
##'        each of dimension d
##' @param signs corresponding matrix of signs sign(x_1), .., sign(x_n)
##' @param l.off the offset to subtract and re-add; ideally in the order of max(.)
##' @param strict logical indicating if it should stop on some negative sums
##' @return log(x_1 + .. + x_n) [i.e., of dimension d] computed via
##'         log(sum(x)) = log(sum(sign(x)*|x|)) = log(sum(sign(x)*exp(log(|x|))))
##'         = log(exp(log(x0))*sum(signs*exp(log(|x|)-log(x0))))
##'         = log(x0) + log(sum(signs* exp(log(|x|)-log(x0))))
##'         = l.off   + log(sum(signs* exp(lxabs -  l.off  )))
##' @author Marius Hofert and Martin Maechler
lssum <- function (lxabs, signs, l.off = apply(lxabs, 2, max), strict = TRUE) {
    stopifnot(length(dim(lxabs)) == 2L) # is.matrix(.) generalized
    sum. <- colSums(signs * exp(lxabs - rep(l.off, each=nrow(lxabs))))
    if (any(is.nan(sum.) || sum. <= 0))
        (if(strict) stop else warning)("lssum found non-positive sums")
    l.off + log(sum.)
}


##' Compute Stirling numbers of the 1st kind
##'
##' s(n,k) = (-1)^{n-k} times
##' the number of permutations of 1,2,...,n with exactly k cycles
##'
##' NIST DLMF 26.8 --> http://dlmf.nist.gov/26.8
##'
##' @title  Stirling Numbers of the 1st Kind
##' @param n
##' @param k
##' @return s(n,k)
##' @author Martin Maechler
Stirling1 <- function(n,k)
{
    ## NOTA BENE: There's no "direct" method available here
    stopifnot(length(n) == 1, length(k) == 1)
    if (k < 0 || n < k) stop("'k' must be in 0..n !")
    if(n == 0) return(1)
    if(k == 0) return(0)
    S1 <- function(n,k) {
        if(k == 0 || n < k) return(0)
        if(is.na(S <- St[[n]][k])) {
            ## s(n,k) = s(n-1,k-1) - (n-1) * s(n-1,k) for all n, k >= 0
            St[[n]][k] <<- S <- if(n1 <- n-1L)
                S1(n1, k-1) - n1* S1(n1, k) else 1
        }
        S
    }
    if(compute <- (nt <- length(St <- get("S1.tab", .nacopEnv))) < n) {
        ## extend the "table":
        length(St) <- n
        for(i in (nt+1L):n) St[[i]] <- rep.int(NA_real_, i)
    }
    else compute <- is.na(S <- St[[n]][k])
    if(compute) {
        S <- S1(n,k)
        ## store it back:
        assign("S1.tab", St, envir = .nacopEnv)
    }
    S
}

##' Full Vector of Stirling Numbers of the 1st Kind
##'
##' @title  Stirling1(n,k) for all k = 1..n
##' @param n
##' @return the same as sapply(1:n, Stirling1, n=n)
##' @author Martin Maechler
Stirling1.all <- function(n)
{
    stopifnot(length(n) == 1)
    if(!n) return(numeric(0))
    if(get("S1.full.n", .nacopEnv) < n) {
        assign("S1.full.n", n, envir = .nacopEnv)
        unlist(lapply(seq_len(n), Stirling1, n=n))# which fills "S1.tab"
    }
    else get("S1.tab", .nacopEnv)[[n]]
}

##' Compute Stirling numbers of the 2nd kind
##'
##' S^{(k)}_n = number of ways of partitioning a set of $n$ elements into $k$
##'	non-empty subsets
##' (Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
##'   Closed Form : p.824 "C."
##'
##' @title  Stirling Numbers of the 2nd Kind
##' @param n
##' @param k
##' @param method
##' @return S(n,k) = S^{(k)}_n
##' @author Martin Maechler, "direct": May 28 1992
Stirling2 <- function(n,k, method = c("lookup.or.store","direct"))
{
    stopifnot(length(n) == 1, length(k) == 1)
    if (k < 0 || n < k) stop("'k' must be in 0..n !")
    method <- match.arg(method)
    switch(method,
           "direct" = {
               sig <- rep(c(1,-1)*(-1)^k, length.out= k+1) # 1 for k=0; -1 1 (k=1)
               k <- 0:k # (!)
               ga <- gamma(k+1)
               round(sum( sig * k^n /(ga * rev(ga))))
           },
           "lookup.or.store" = {
               if(n == 0) return(1) ## else:
               if(k == 0) return(0)
               S2 <- function(n,k) {
                   if(k == 0 || n < k) return(0)
                   if(is.na(S <- St[[n]][k]))
                       ## S(n,k) = S(n-1,k-1) + k * S(n-1,k)   for all n, k >= 0
                       St[[n]][k] <<- S <- if(n1 <- n-1L)
                           S2(n1, k-1) + k* S2(n1, k) else 1 ## n = k = 1
                   S
               }
               if(compute <- (nt <- length(St <- get("S2.tab", .nacopEnv))) < n) {
                   ## extend the "table":
                   length(St) <- n
                   for(i in (nt+1L):n) St[[i]] <- rep.int(NA_real_, i)
               }
               else compute <- is.na(S <- St[[n]][k])
               if(compute) {
                   S <- S2(n,k)
                   ## store it back:
                   assign("S2.tab", St, envir = .nacopEnv)
               }
               S
           })
}

##' Full Vector of Stirling Numbers of the 2nd Kind
##'
##' @title  Stirling2(n,k) for all k = 1..n
##' @param n
##' @return the same as sapply(1:n, Stirling2, n=n)
##' @author Martin Maechler
Stirling2.all <- function(n)
{
    stopifnot(length(n) == 1)
    if(!n) return(numeric(0))
    if(get("S2.full.n", .nacopEnv) < n) {
        assign("S2.full.n", n, envir = .nacopEnv)
        unlist(lapply(seq_len(n), Stirling2, n=n))# which fills "S2.tab"
    }
    else get("S2.tab", .nacopEnv)[[n]]
}


##' Compute Eulerian numbers (German "Euler Zahlen")  A(n,k)
##'
##' A(n,k) = number of permutations of n with exactly  k ascents (or k descents)
##' --> http://dlmf.nist.gov/26.14
##'
##' @title Eulerian Numbers
##' @param n
##' @param k
##' @param method
##' @return A(n,k) = < n \\ k >
##' @author Martin Maechler, April 2011
Eulerian <- function(n,k, method = c("lookup.or.store","direct"))
{
    stopifnot(length(n) == 1, length(k) == 1)
    if(k < 0 || n < k) stop("'k' must be in 0..n !")
    if(n && k == n) return(0)
    ## have  __ 0 <= k < n __
    method <- match.arg(method)
    switch(method,
	   "direct" = { ## suffers from cancellation eventually ..
	       if(k == 0) return(1)
	       if(k == n) return(0)
	       ## else 0 <= k < n >= 2
	       ## http://dlmf.nist.gov/26.14.E9 :  A(n,k) = A(n, n-1-k),  n >= 1
	       if(k >= (n+1)%/% 2) k <- n-(k+1L)
	       k1 <- k+1L
	       sig <- rep(c(1,-1), length.out = k1) # 1 for k=0; 1 -1 (k=1), ...
	       round(sum( sig * choose(n+1, 0:k) * (k1:1L)^n ))
	   },
	   "lookup.or.store" = {
	       Eul <- function(n,k) {
		   ## Quick return for those that are *not* stored:
		   if(k < 0 || k > n || (0 < n && n == k)) return(0)
		   if(n == 0) return(1)
		   ## now -- 0 <= k < n -- are stored
		   if(is.na(r <- E.[[n]][k1 <- k+1L])) ## compute it (via recursion)
		       ## A(n,k) = (k+1)* A(n-1,k) + (n-k)*A(n-1,k-1)	for n >= 2, k >= 0
		       ## http://dlmf.nist.gov/26.14.E8
		       E.[[n]][k1] <<- r <- if(n1 <- n-1L)
			   k1*Eul(n1, k)+ (n-k)* Eul(n1, k-1) else 1 ## n=1, k=0
		   r
	       }
	       if(compute <- (nt <- length(E. <- get("Eul.tab", .nacopEnv))) < n) {
		   ## extend the "table":
		   length(E.) <- n
		   for(i in (nt+1L):n) E.[[i]] <- rep.int(NA_real_, i)
	       }
	       else compute <- is.na(E <- E.[[n]][k+1L])
	       if(compute) {
		   E <- Eul(n,k)
		   ## store it back:
		   assign("Eul.tab", E., envir = .nacopEnv)
	       }
	       E
	   })
}

##' Full Vector of Eulerian Numbers == row of Euler triangle
##'
##' @title Eulerian(n,k) for all k = 0..n-1
##' @param n
##' @return (for n >= 1), the same as sapply(0:(n-1), Eulerian, n=n)
##' @author Martin Maechler, April 2011
Eulerian.all <- function(n)
{
    stopifnot(length(n) == 1, n >= 0)
    if(!n) return(1)
    if(get("Eul.full.n", .nacopEnv) < n) {
        ## FIXME: do the assign() only when the lapply() below does *not* fail
        ##  on.exit( ... ) ?
	assign("Eul.full.n", n, envir = .nacopEnv)
	unlist(lapply(0:(n-1L), Eulerian, n=n))# which fills "Eul.tab"
    }
    else get("Eul.tab", .nacopEnv)[[n]]
}

## Our environment for tables etc:  no hash, as it will contain *few* objects:
.nacopEnv <- new.env(parent=emptyenv(), hash=FALSE)
assign("S2.tab", list(), envir = .nacopEnv) ## S2.tab[[n]][k] == S(n, k)
assign("S1.tab", list(), envir = .nacopEnv) ## S1.tab[[n]][k] == s(n, k)
assign("S2.full.n", 0  , envir = .nacopEnv)
assign("S1.full.n", 0  , envir = .nacopEnv)
assign("Eul.tab", list(), envir = .nacopEnv) ## Eul.tab[[n]][k] == A(n, k) == < n \\ k > (Eulerian)
assign("Eul.full.n", 0	, envir = .nacopEnv)


##' From   http://en.wikipedia.org/wiki/Polylogarithm
##' 1. For integer values of the polylogarithm order, the following
##'   explicit expressions are obtained by repeated application of  z * d/dz
##'   to Li_1(z):
##' ---
##'     {Li}_{1}(z) = -\ln(1-z)
##'     {Li}_{0}(z) = {z \over 1-z}
##'     {Li}_{-1}(z) = {z \over (1-z)^2}
##'     {Li}_{-2}(z) = {z \,(1+z) \over (1-z)^3}
##'     {Li}_{-3}(z) = {z \,(1+4z+z^2) \over (1-z)^4}
##'     {Li}_{-4}(z) = {z \,(1+z) (1+10z+z^2) \over (1-z)^5}.
##' ---
##' Accordingly the polylogarithm reduces to a ratio of polynomials in
##' z, and is therefore a rational function of z, for all nonpositive
##' integer orders. The general case may be expressed as a finite sum:
##' ---
##' {Li}_{-n}(z) = \left(z \,{\partial \over \partial z} \right)^n \frac{z}{1-z} =
##'     = \sum_{k=0}^n k! \,S(n+1,k+1) \left({z \over {1-z}} \right)^{k+1}
##' \ \ (n=0,1,2,\ldots),
##' ---
##' where S(n,k) are the Stirling numbers of the second kind.
##' Equivalent formulae applicable to negative integer orders are
##' (Wood 1992, § 6):
##' ---
##'  {Li}_{-n}(z) = (-1)^{n+1} \sum_{k=0}^n k! S(n+1,k+1) (\frac{-1}{1-z})^{k+1} \
##'     (n=1,2,3,...),
##' and:
##'  {Li}_{-n}(z) = \frac{1}{(1-z)^{n+1}} sum_{k=0}^{n-1} < n \ k >  z^{n-k}
##'               = \frac{z \sum_{k=0}^{n-1} < n \ k >  z^k} {(1-z)^{n+1}},
##' where  < n \ k >  are the  Eulerian numbers.
##' All roots of Li_n(z) are distinct and real; they include z = 0.
##' Duplication formula:  2^{1-s} Li_s(z^2) = Li_s(z) + Li_s(-z).
##'
##' Compute the polylogarithm function \eqn{Li_s(z)}
##'
##' @title Polylogarithm Li_s(z)
##' @param z numeric or complex vector
##' @param s complex number; current implementation is aimed at s \in 0,-1,..
##' @param method a string specifying the algorithm to be used
##' @param logarithm logical specifying to return log(Li.(.)) instead of Li.(.)
##' @param is.log.z logical; if TRUE, the specified 'z' argument is really  w = log(z);
##' 	i.e., compute  Li_s(exp(w))  {and we typically have w < 0  <==> z < 1}
##' @param n.sum  for "sum"--this is more for experiments etc
##' @return numeric/complex vector as \code{z}
##' @author Martin Maechler
polylog <- function(z, s, method = c("default", "sum", "negI-s-Stirling",
			  "negI-s-Eulerian", "negI-s-asymp-w"),
		    logarithm = FALSE, is.log.z = FALSE,
                    is.logmlog = FALSE, asymp.w.order = 0,
		    ## for "sum" -- this is more for experiments etc:
		    n.sum)
{
    if((nz <- length(z)) == 0 || (ns <- length(s)) == 0)
	return((z+s)[FALSE])# of length 0
    stopifnot(length(s) == 1) # for now
    method <- match.arg(method)
    if(is.log.. <- is.log.z) {
	if(is.logmlog)
	    stop("cannot have *both* 'is.log.z' and 'is.logmlog'")
	## z = e^{w}
	w <- z
    } else if(is.log.. <- is.logmlog) {
	stopifnot((lw <- z) < 0)
	w <- -exp(lw)
    }
    if(is.log..) z <- exp(w)

    if(s == 2 && method == "default")## Dilog aka Dilogarithm, from gsl pkg -> GSL = GNU Scientific Library
        return(if(is.numeric(z)) dilog(z) else complex_dilog(z))
    if(method == "default")
        method <- if(s == as.integer(s) && s <= 1)
                      "negI-s-Stirling" else "sum"
    if(method == "sum") {
	stopifnot((Mz <- Mod(z)) <= 1, Mz < 1 | Re(s) > 1,
		  n.sum > 99, length(n.sum) == 1)
	p <- polynEval((1:n.sum)^-s, z)
	if(logarithm) (if(is.log..) w else log(z)) +log(p) else z*p
    } else {
	## "negI"-methods:  s \in {1, 0, -1, -2, ....}	 "[neg]ative [I]nteger"
	stopifnot(s == as.integer(s), s <= 1)
	if(s == 1) { ## result = -log(1 -z) = -log(1 - exp(w)) = -log1mexp(-w)
	    r <- if(is.log..) -log1mexp(-w) else -log1p(-z)
	    return(if(logarithm) log(r) else r)
	}
	## s <= 0:

	## i.z := (1 - z) = 1 - exp(w)
	i.z <- if(is.log.z) -expm1(w) else (1 - z)
	n <- -as.integer(s)
	switch(method,
	       "negI-s-Stirling" =
	   {
	       r <- z/i.z ## = z/(1 - z)
	       ## if(n == 0) return(r)
	       ## k1 <- seq_len(n+1)# == k+1, k = 0...n
	       fac.k <- cumprod(c(1, seq_len(n)))
	       S.n1.k1 <- Stirling2.all(n+1) ## == Stirling2(n+1, k+1)
	       p <- polynEval(fac.k * S.n1.k1, r)
	       ## log(r) = log(z / (1-z)) = logit(z) = qlogis(z),  and if(is.log.z)
	       ## log(r) = log(z / (1-z)) = log(z) - log(1-z) = w - log(1-exp(w)) =
	       ##	 = w - log1mexp(-w)
	       if(logarithm) log(p) + if(is.log..) w - log1mexp(-w) else qlogis(z)
	       else r*p
	   },
	       "negI-s-Eulerian" =
	   {
	       ## Li_n(z) = \frac{z \sum_{k=0}^{n-1} < n \ k >	z^k} {(1-z)^{n+1}},
	       Eu.n <- Eulerian.all(n)
	       ## FIXME? do better for z=exp(w), or large n --> ## log(Eulerian) etc
	       p <- polynEval(Eu.n, z)
	       if(logarithm)
		   log(p) +
		       if(is.log..) w - (n+1)*log1mexp(-w)
		       else log(z) - (n+1)*log1p(-z)
	       else z*p / i.z^(n+1)
	   },
	       "negI-s-asymp-w" =
	       ## MM's asymptotic formula for small w = log(z):
	   {
	       stopifnot(length(asymp.w.order) == 1, asymp.w.order >= 0)
	       if(!is.numeric(z) || z < 0)
		   stop("need z > 0, or rather w = log(z) < 0 for method ",
			method,"; but z=",z)
	       if(!is.log..) w <- log(z)
	       w <- -w # --> w > 0 , z = exp(-w)  as in "paper", lw = log(w), w = exp(lw)
	       if(w >= 1)
		   warning("w << 1 is prerequisite for method ",
			   method, "; but w=", w)
	       n1 <- n+1L
	       if(asymp.w.order == 0) { ## zero-th-order ("simple")
		   if(logarithm) lgamma(n1) - n1*(if(is.logmlog)lw else log(w))
		   else gamma(n1)* if(is.logmlog) exp(-lw*n1) else w^-n1
	       } else if(asymp.w.order == 1) { ## 1st order
		   stop("1-st order asympt: implementation unfinished")
		   r <- if(logarithm) {
		       ## main: lgamma(n+1) - (n+1)*log(w) + c1 ^ w^*(2*ceiling((n+1)/2))
		       trm <-
			   if(n %% 2) { ## n odd
			   } else { ## n even
			   }
		   } else gamma(n+1)*w^(-(n+1)) +0-0+0-0+ Bernoulli()/factorial()

	       } else stop("asymp.w.order = ",asymp.w.order," is not implemented")

### TODO:
###  Frank: large theta: even w underflows to 0, ---> use  is.logmlog = TRUE
###
### TODO 2): use *first* order term  ... [ for small w  i.e. z ~<~ 1 ]
	   },

	       stop("unsupported method ", method))
    }
}


### FIXME: Do *cache* the full triangle,  "similarly" to Stirling
##  ------

##' Generate all Bernoulli numbers up to the n-th,
##' using diverse algorithms
##'
##' @title Bernoulli Numbers up to Order n
##' @param n
##' @param method
##' @param precBits
##' @param verbose
##' @return numeric vector of length n, containing B(n)
##' @author Martin Maechler, 25 Jun 2011 (night train Vienna - Zurich)
Bernoulli.all <-
    function(n, method = c("A-T", "sumBin", "sumRamanujan", "asymptotic"),
             precBits = NULL, verbose = getOption("verbose"))
{
    stopifnot(length(n <- as.integer(n)) == 1)
    method <- match.arg(method)
    switch(method,
           "A-T" =
       {
           if(verbose) stopifnot(requireNamespace("MASS"))
           nn <- seq_len(n1 <- n+1L) ## <- FIXME, make this work with Rmpfr optionally
           if(!is.numeric(precBits) || !is.finite(precBits)) {
               B <- numeric(n1)
               Bv <- 1/nn
           } else {
	       stopifnot(requireNamespace("Rmpfr", quietly=TRUE),
			 length(precBits) == 1, precBits >= 10)
	       mpfr <- Rmpfr::mpfr
               B <- mpfr(numeric(n1), precBits=precBits)
               Bv <- mpfr(1, precBits=precBits)/nn
           }
           for(m in seq_len(n1)) {
               ## length(Bv) == n+1+1-m = n+2-m
               if(verbose) { cat(m-1,":"); print(MASS::fractions(Bv)) }
               B[m] <- Bv[1L]
               if(m <= n) ## recursion:
                   Bv <- seq_len(n1-m) * (Bv[-(n1+1L-m)] - Bv[-1L])
           }
       },
           "sumBin" = , "sumRamanujan" = , "asymptotic" =
       {
           B <- c(1, -1/2)
           if(n <= 1) return(B[1:(n+1)])
           Bn <- Bernoulli(n, method=method)#-> fill "Bern.tab" cache
           B <- c(B, numeric(n-2), Bn)
           ##     0:1  2...n-1
           ii <- seq_len((n-1L) %/% 2)
           B[ii*2L] <- get("Bern.tab", .nacopEnv)[[method]][ii]
       },
           stop(gettextf("method '%s' not yet implemented", method), domain=NA)
           )## end{switch}
    B
}

## FIXME: The Akiyama-Tanigawa algorithm formula is numerically >> bad <<

## Instead, use the recursive definition of Bernoulli numbers, or

assign("Bern.tab", list(), envir = .nacopEnv) ## Bern.tab[[method]][n] == Bernoulli_{2n}

Bernoulli <- function(n, method = c("sumBin", "sumRamanujan", "asymptotic"),
                      verbose = FALSE)
{
    if(n <= 1) {
        if(n == 0) 1 else -1/2
    } else if(n %% 2) 0
    else { ## n >= 2, even
	method <- match.arg(method)
	## if(double.precision)
	if(n > 224) return(if(n %% 4) -Inf else Inf) # overflow
	n2 <- n %/% 2
	if(length(Bt <- get("Bern.tab", .nacopEnv)[[method]]) < n2 ||
	   is.na(B <- Bt[n2])) {
	    ## compute it -- and cache it
	    switch(method,
                   "asymptotic" =
               {
                   ## B <- (-1)^(n2+1) * 2*(1 + 2^-n)* exp(lfactorial(n) - n * log(2*pi))
                   B <- (-1)^(n2+1) * 2*exp(lgamma(n+1) + log1p(2^-n + 3^-n) - n * log(2*pi))
                   ## where we replace zeta(n) = 1 + 1/2^n + 1/3^n + 1/4^n + ...
                   ## by   1 + 1/2^n + 1/3^n
                   ## http://en.wikipedia.org/wiki/Bernoulli_number#Asymptotic_approximation
               },
		   "sumBin" =
	       {
		   ## This is R-slow for n ~= 1000, but ok for n ~= 100
		   s <- 0
		   binc <- 1
		   n.1 <- n - 1L
		   for(j in 0:n.1) { ## binc == choose(n, j)
                       if(verbose) {
                           cat(sprintf("(n,j)=(%d,%d): Bin.cf= %g", n,j, round(binc)))
                           Bj <- Bernoulli(j)
                           cat(sprintf(" Bj= %g\n", Bj))
                           s <- s + binc * Bj / (n-j+1L)
                       } else s <- s + binc * Bernoulli(j) / (n-j+1L)
		       if(j < n.1)
			   binc <- round(binc * (n-j)/(j+1L))
		   }
		   B <- -s
		   ## Cache it: update the *current* table (updated in the mean time!)
		   Btab <- get("Bern.tab", .nacopEnv)
		   Btab[[method]][n2] <- B
		   assign("Bern.tab", Btab, envir = .nacopEnv)
	       },
		   stop(gettextf("method '%s' not yet implemented", method), domain=NA)
		   )## end{switch}
	}
	B
    }
}


## slightly more efficiently,  Ramanujan's congruences:
## == http://en.wikipedia.org/wiki/Bernoulli_number#Ramanujan.27s_congruences
##------------------------------------------------------------------------
## Ramanujan's congruences

## The following relations, due to Ramanujan, provide a more efficient method for calculating Bernoulli numbers:

##     {{m+3}\choose{m}}B_m=\begin{cases} {{m+3}\over3}-\sum\limits_{j=1}^{m/6}{m+3\choose{m-6j}}B_{m-6j}, & \mbox{if}\ m\equiv 0\pmod{6};\\ {{m+3}\over3}-\sum\limits_{j=1}^{(m-2)/6}{m+3\choose{m-6j}}B_{m-6j}, & \mbox{if}\ m\equiv 2\pmod{6};\\ -{{m+3}\over6}-\sum\limits_{j=1}^{(m-4)/6}{m+3\choose{m-6j}}B_{m-6j}, & \mbox{if}\ m\equiv 4\pmod{6}.\end{cases}
##------------------------------------------------------------------------

##---> Free Python implementations are here
##>> http://en.literateprograms.org/Bernoulli_numbers_%28Python%29
##--> MM: see ../misc/Bernoulli_numbers_(Python)/


debye1 <- function(x) {
    ## gsl's debye_1() gives *errors* on NaN/NA/Inf
    if(any(nfin <- !(fin <- is.finite(x)))) {
	d <- abs(x)
	d[fin] <- debye_1(d[fin])
	if(any(isI <- d[nfin] == Inf))
	    d[nfin][which(isI)] <- 0 # which(.) drops NAs
	if(length(neg <- which(x < 0))) # NA's !
	    d[neg] <- d[neg] - x[neg]/2
	d
    }
    else {
	d <- debye_1(abs(x))
	## ifelse(x >= 0, d, d - x / 2) ## k = 1, Frees & Valdez 1998, p.9
	d - (x<0) * x / 2
    }
}


debye2 <- function(x) {
    ## gsl's debye_2() gives *errors* on NaN/NA/Inf
    if(any(nfin <- !(fin <- is.finite(x)))) {
	d <- abs(x)
	d[fin] <- debye_2(d[fin])
	if(any(isI <- d[nfin] == Inf))
	    d[nfin][which(isI)] <- 0 # which(.) drops NAs
	if(length(neg <- which(x < 0))) # NA's !
	    d[neg] <- d[neg] - 2/3 * x[neg]
	d
    }
    else {
	## k = 2, Frees & Valdez 1998, p.9
	d <- debye_2(abs(x))
	d - 2/3 * (x<0) * x
    }
}
