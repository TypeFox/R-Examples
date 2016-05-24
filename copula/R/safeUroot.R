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


##' @title (Safely) Finding a root
##' @param f function
##' @param interval interval
##' @param ...
##' @param lower lower endpoint
##' @param upper upper endpoint
##' @param f.lower function value at lower endpoint
##' @param f.upper function value at upper endpoint
##' @param Sig *desired* sign of f(upper), or NULL
##' @param tol tolerance
##' @param maxiter maximal number of iterations
##' @param trace number determining tracing
##' @author Martin Maechler (from Martin's package "nor1mix")
safeUroot <-
    function(f, interval, ...,
	     lower = min(interval), upper = max(interval),
	     f.lower = f(lower, ...), f.upper = f(upper, ...),
	     Sig = NULL, check.conv = FALSE,
	     tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
{
    if(	  is.null(Sig) && f.lower * f.upper > 0 ||
       is.numeric(Sig) && (Sig*f.lower > 0 || Sig*f.upper < 0)) {
	if(trace)
	    cat(sprintf("search in [%g,%g]%s", lower, upper,
			if(trace >= 2)"\n" else " ... "))
	Delta <- function(u) 0.01* pmax(1e-7, abs(u))
	## Two cases:
	if(is.null(Sig)) {
	    ## case 1)	'Sig' unspecified --> extend (lower, upper) at the same time
	    delta <- Delta(c(lower,upper))
	    while(isTRUE(f.lower*f.upper > 0) && any(iF <- is.finite(c(lower,upper)))) {
		if(iF[1]) f.lower <- f(lower <- lower - delta[1])
		if(iF[2]) f.upper <- f(upper <- upper + delta[2])
		if(trace >= 2)
		    cat(sprintf(" .. modified lower,upper: (%15g,%15g)\n",
				lower,upper))
		delta <- 2 * delta
	    }
	} else {
	    ## case 2) 'Sig' specified --> typically change only *one* of lower, upper
	    ## make sure we have Sig*f(lower) < 0 and Sig*f(upper) > 0:
	    delta <- Delta(lower)
	    while(isTRUE(Sig*f.lower > 0)) {
		f.lower <- f(lower <- lower - delta)
		if(trace) cat(sprintf(" .. modified lower: %g\n",lower))
		delta <- 2 * delta
	    }
	    delta <- Delta(upper)
	    while(isTRUE(Sig*f.upper < 0)) {
		f.upper <- f(upper <- upper + delta)
		if(trace) cat(sprintf(" .. modified upper: %g\n",upper))
		delta <- 2 * delta
	    }
	}
	if(trace && trace < 2)
	    cat(sprintf("extended to [%g, %g]\n", lower, upper))
    }
    if(!isTRUE(f.lower * f.upper <= 0))
	stop("did not succeed extending the interval endpoints for f(lower) * f(upper) <= 0")
    if(check.conv) {
	r <- tryCatch(uniroot(f, ..., lower=lower, upper=upper,
			      f.lower = f.lower, f.upper = f.upper,
			      tol=tol, maxiter=maxiter),
		      warning = function(w)w)
	if(inherits(r,"warning"))
	    stop("convergence problem in zero finding: ", r$message)
	else r
    }
    else
	uniroot(f, ..., lower=lower, upper=upper,
		f.lower = f.lower, f.upper = f.upper,
		tol=tol, maxiter=maxiter)
}
