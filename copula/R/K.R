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


### Kendall distribution #######################################################

## deprecated (former) Kendall distribution function K
K <- function(u, cop, d, n.MC=0, log=FALSE){
    .Deprecated("pK") # set K to "deprecated" => throws a message
    pK(u, cop=cop, d=d, n.MC=n.MC, log.p=log) # call the new function
}

##' Empirical Kendall distribution function K_{n,d} as in Lemma 1 of
##' Genest, Neslehova, Ziegel (2011)
##'
##' @title Empirical Kendall distribution function
##' @param u evaluation points u in [0,1]
##' @param x data (in IR^d) based on which K is estimated
##' @return K_{n,d}(u)
##' @author Marius Hofert
##' Note: This is the empirical distribution function of a discrete radial part
##'       and thus K_{n,d}(0) > 0. The mass at the largest value of the support
##'       of R determines K_{n,d}(0)
Kn <- function(u, x)
{
    stopifnot(0 <= u, u <= 1, (n <- nrow(x)) >= 1, (d <- ncol(x)) >= 1)
    W <- vapply(seq_len(n), function(i) sum( colSums(t(x)<x[i,])==d ) / (n+1), NA_real_)
    ecdf(W)(u)
}

##' Distribution function of the Kendall distribution
##'
##' @title Kendall distribution function
##' @param u evaluation point(s) in [0,1]
##' @param cop acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal
##'	   to n.MC to evaluate the generator derivatives; otherwise the exact
##'        formula is used
##' @param log.p logical indicating whether the logarithm is returned
##' @return Kendall distribution function at u
##' @author Marius Hofert
pK <- function(u, cop, d, n.MC=0, log.p=FALSE)
{
    stopifnot(is(cop, "acopula"), 0 <= u, u <= 1)
    ## limiting cases
    n <- length(u)
    res <- numeric(n)
    res[is0 <- u == 0] <- if(log.p) -Inf else 0
    res[is1 <- u == 1] <- if(log.p) 0 else 1
    not01 <- seq_len(n)[!(is0 | is1)]
    uN01 <- u[not01]
    ## computations
    th <- cop@theta
    if(n.MC > 0) { # Monte Carlo
	stopifnot(is.finite(n.MC))
        if(length(not01)){
            V <- cop@V0(n.MC, th) # vector of length n.MC
            psiI <- cop@iPsi(uN01, th) # vector of length n
            lr <- unlist(lapply(psiI, function(psInv){
                -log(n.MC) + lsum(as.matrix(ppois(d-1, V*psInv, log.p=TRUE)))
                ## Former code: mean(ppois(d-1, V*psInv))
            }))
            res[not01] <- if(log.p) lr else exp(lr)
        }
    } else { # direct
	if(length(not01))
	    res[not01] <- if(d == 1) { # d == 1
		if(log.p) log(uN01) else uN01 # K(u) = u
	    } else if(d == 2) { # d == 2
		r <- uN01 + exp( cop@iPsi(uN01, theta=th, log=TRUE) -
                                 cop@absdiPsi(uN01, th, log=TRUE) ) # K(u) = u - psi^{-1}(u) / (psi^{-1})'(u)
                if(log.p) log(r) else r
	    } else { # d >= 3
		j <- seq_len(d-1)
		lpsiI. <- cop@iPsi(uN01, theta=th, log=TRUE)
		labsdPsi <- do.call(rbind,
				    lapply(j, function(j.)
					   cop@absdPsi(exp(lpsiI.),
						       theta=th,
						       degree=j.,
						       log=TRUE))) # (d-1) x n  matrix [n = length(not01)]
                ## containing log( (-1)^j * psi^{(j)}(psi^{-1}(u)) ) in the j-th row
		lfac.j <- cumsum(log(j)) ## == lfactorial(j)
                lx <- labsdPsi + j %*% t(lpsiI.) - lfac.j # (d-1) x n matrix
                lx <- rbind(log(uN01), lx) # d x n matrix containing the logarithms of the summands of K
                ls <- lsum(lx) # log(K(u))
                if(log.p) ls else pmin(1, exp(ls)) # ensure we are in [0,1] {numerical inaccuracy}
		## Former code:
		## K2 <- function(psInv) {
		##    labsdPsi <- unlist(lapply(j, cop@absdPsi,
		##			      u=psInv, theta=th, log=TRUE))
		##    sum(exp(labsdPsi + j*log(psInv) - lfac.j))
		## }
		## pmin(1, uN01 + unlist(lapply(psiI[not01], K2)))
		##
		## NB: AMH, Clayton, Frank are numerically not quite monotone near one;
		## --  this does not change that {but maybe slightly *more* accurate}:
		## absdPsi. <- unlist(lapply(j, cop@absdPsi, u = psInv, theta = th,
		##						 log = FALSE))
		##		       sum(absdPsi.*psInv^j/factorial(j))
	    } # else (d >= 3)
    } # if/else method (MC/direct)
    res
}


##' Quantile function of the Kendall distribution
##'
##' @title Quantile function of the Kendall distribution function
##' @param u evaluation point(s) in [0,1]
##' @param cop acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied to evaluate K with
##'        sample size equal to n.MC; otherwise the exact formula is used
##' @param method method used for inverting K; currently:
##'        "default" : chooses a useful default
##'        "simple"  : straightforward root finding
##'        "sort"    : root finding after sorting the u's
##'        "discrete": evaluating K on u.grid, then finding approximating
##'                    quantiles based on these values
##'        "monoH.FC": evaluating K on u.grid, then approximating K via monotone
##'                    splines; finding quantiles via uniroot
##' @param u.grid default grid on which K is computed for method "discrete" and
##'        "monoH.FC"
##' @param ... additional arguments passed to uniroot() (for methods "sort",
##'        "simple", and "monoH.FC") or findInterval() (for method "discrete")
##' @return Quantile function of the Kendall distribution function at u
##' @author Marius Hofert
##' Note: - K(u) >= u => u gives an upper bound for K^{-1}(u)
##'       - K for smaller dimensions would also give upper bounds for K^{-1}(u),
##'         but even for d=2, there is no explicit formula for K^{-1}(u) known.
qK <- function(u, cop, d, n.MC=0,
               method=c("default", "simple", "sort", "discrete", "monoH.FC"),
               u.grid, ...)
{
    stopifnot(is(cop, "acopula"), 0 <= u, u <= 1)
    if(d==1) ## special case : K = identity
	return(u)
    ## limiting cases
    n <- length(u)
    res <- numeric(n)                   # all 0
    res[is1 <- u == 1] <- 1
    is0 <- u == 0 ## res[is0 <- u == 0] <- 0
    if(!any(not01 <- !(is0 | is1)))
	return(res)
    ## usual case:
    uN01 <- u[not01]           # u's for which we have to determine quantiles
    lnot01 <- sum(not01)       # may be < n
    ## computing the quantile function
    method <- match.arg(method)
    res[not01] <-
        switch(method,
               "default" =
           {
               ## Note: This is the same code as method="monoH.FC" (but with a
               ##       chosen grid)
               u.grid <- 0:128/128 # default grid
               K.u.grid <- pK(u.grid, cop=cop, d=d, n.MC=n.MC, log.p=FALSE)
               ## function for root finding
               fspl <- function(x, u)
                   splinefun(u.grid, K.u.grid, method = "monoH.FC")(x) - u
               ## root finding
               vapply(uN01, function(u)
                      uniroot(fspl, u=u, interval=c(0,u), ...)$root, NA_real_)

           },
               "simple" =               # straightforward root finding
           {
               ## function for root finding
               f <- function(t, u) pK(t, cop=cop, d=d, n.MC=n.MC, log.p=FALSE) - u
               ## root finding
               vapply(uN01, function(u)
                      uniroot(f, u=u, interval=c(0,u), ...)$root, NA_real_)
           },
               "sort" =                 # root finding with first sorting u's
           {
               ## function for root finding
               f <- function(t, u) pK(t, cop=cop, d=d, n.MC=n.MC, log.p=FALSE) - u
               ## sort u's
               ord <- order(uN01, decreasing=TRUE)
               uN01o <- uN01[ord]       # uN01 ordered in decreasing order
               ## deal with the first one (largest u) separately
               q <- numeric(lnot01)
               q[1] <- uniroot(f, u=uN01o[1], interval=c(0, uN01o[1]), ...)$root # last quantile -- used in the following
               if(lnot01 > 1){
                   for(i in 2:lnot01){
                       lower <- 0
                       eps <- 1e-4 # ugly but due to non-monotonicity of K [otherwise: "Error... f() values at end points not of opposite sign"]
                       upper <- min(uN01o[i], q[i-1]+eps)
                       if(FALSE){       # checks for debugging
                           f.lower <- f(lower, uN01o[i])
                           f.upper <- f(upper, uN01o[i])
                           if(lower >= upper) stop("lower=", lower, ", upper=", upper)
                           if(sign(f.lower*f.upper) >= 0) stop("uN01o[",i,"]=", uN01o[i], ", f.lower=", f.lower, ", f.upper=", f.upper)
                       }
                       q[i] <- uniroot(f, u=uN01o[i], interval=c(lower, upper), ...)$root
                   }
               }
               ## "return" unordered result
               quo <- numeric(lnot01)   # unordered vector of quantiles q
               quo[ord] <- q            # return unordered result
               quo
           },
               "discrete" = # evaluate K at a grid and compute approximate quantiles based on this grid
           {
               stopifnot(0 <= u.grid, u.grid <= 1)
               K.u.grid <- pK(u.grid, cop=cop, d=d, n.MC=n.MC, log.p=FALSE)
               u.grid[findInterval(uN01, vec=K.u.grid,
                                   rightmost.closed=TRUE, ...)+1] # note: this gives quantiles according to the "typical" definition
           },
               "monoH.FC" = # root finding based on an approximation of K via monotone splines (see ?splinefun)
           {
               ## evaluate K at a grid
               stopifnot(0 <= u.grid, u.grid <= 1)
               K.u.grid <- pK(u.grid, cop=cop, d=d, n.MC=n.MC,
                              log.p=FALSE)
               ## function for root finding
               fspl <- function(x, u)
                   splinefun(u.grid, K.u.grid, method = "monoH.FC")(x) - u
               ## root finding
               vapply(uN01, function(u)
                      uniroot(fspl, u=u, interval=c(0,u), ...)$root, NA_real_)
           },
               stop("unsupported method ", method))
    res
}


##' Density of the Kendall distribution
##'
##' @title Density of the Kendall distribution
##' @param u evaluation point(s) in (0,1)
##' @param cop acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal
##'	   to n.MC to evaluate the d-th generator derivative; otherwise the exact
##'        formula is used
##' @param log.p logical indicating whether the logarithm of the density is returned
##' @return Density of the Kendall distribution at u
##' @author Marius Hofert
dK <- function(u, cop, d, n.MC=0, log.p=FALSE)
{
    stopifnot(is(cop, "acopula"), 0 < u, u < 1)
    th <- cop@theta
    lpsiI <- cop@iPsi(u, theta=th, log=TRUE) # log(psi^{-1}(u))
    lpsiIDabs <- cop@absdiPsi(u, theta=th, log=TRUE) # (-psi^{-1})'(u)
    ld <- lfactorial(d-1) # log((d-1)!)
    psiI <- cop@iPsi(u, theta=th) # psi^{-1}(u)
    labsdPsi <- cop@absdPsi(psiI, theta=th, degree=d, n.MC=n.MC, log=TRUE) # log((-1)^d psi^{(d)}(psi^{-1}(u)))
    res <- labsdPsi-ld+(d-1)*lpsiI+lpsiIDabs
    if(log.p) res else exp(res)
}


##' Random number generation for the Kendall distribution
##'
##' @title Random number generation for the Kendall distribution
##' @param n number of random variates to generate
##' @param cop acopula with specified parameter or onacopula
##' @param d dimension (only needed if ...)
##' @return Random numbers from the Kendall distribution
##' @author Marius Hofert
rK <- function(n, cop, d) {
    if(is(cop, c1 <- "acopula")) {
	stopifnot(d == round(d))
	cop <- onacopulaL(cop@name, list(cop@theta, 1L:d))
    } else if(!is(cop, c2 <- "outer_nacopula"))
	stop(gettextf("'cop' must be \"%s\" or \"%s\"", c1,c2), domain=NA)

    pCopula(rCopula(n, cop), cop)
}
