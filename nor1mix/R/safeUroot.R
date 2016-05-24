###  safeUroot()  used to be  *internal*  to qnorMix()
###  -----------  but as it's generally useful,
##   {pay a (very?) small penalty and make it public __ TODO: currently NAMESPACE-hidden!
safeUroot <- function (f, interval, lower = min(interval), upper = max(interval),
		       f.lower = f(lower), f.upper = f(upper),
		       Sig = sign(f.upper),
		       tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0, ...)
{
    if(trace >= 2)
	cat(sprintf("search in [%g,%g]\n", lower, upper))

    ## make sure we have Sig*f(lower) < 0 and Sig*f(upper) > 0:
    delta.r <- 0.01*max(1e-7, abs(lower))
    f.lo <- f.lower
    while(Sig*f.lo > 0) {
	f.lo <- f(lower <- lower - delta.r)
	if(trace)
	    cat(sprintf(" .. modified lower: %g\n",lower))
	delta.r <- 2 * delta.r
    }
    delta.r <- 0.01*max(1e-7, abs(upper))
    f.up <- f.upper
    while(Sig*f.up < 0) {
	f.up <- f(upper <- upper + delta.r)
	if(trace)
	    cat(sprintf(" .. modified upper: %g\n",upper))
	delta.r <- 2 * delta.r
    }

    ## Here, we require R >= 2.6.0 with the improved uniroot():
    uniroot(f, lower=lower, upper=upper,
	    f.lower = f.lo, f.upper = f.up,
	    tol=tol, maxiter=maxiter, ...)
}
