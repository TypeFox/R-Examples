findstysols <-
function(Nsims=100, Ncoefs=3, tsx, tsy, sf=100, plot.it=FALSE,
	print.it=FALSE, verbose=FALSE,
	lctsfn=LCTS, prodcomb.fn=prodcomb, filter.number=1,
	family=c("DaubExPhase", "DaubLeAsymm"), my.maxit=500,
	spec.filter.number=1,
	spec.family=c("DaubExPhase", "DaubLeAsymm"),
	optim.control=list(maxit=my.maxit, reltol=1e-6), irng=rnorm,
	lapplyfn=lapply,Bsims=200, ...)
{
family <- match.arg(family)
spec.family <- match.arg(spec.family)
#
# Find stationary time-varying linear combinations of two time series.
# Produce Nsims solutions, we do this by putting the info into lists
# and use a list processor (eg lapply or mclapply) to perform the calcs
#

#
# Check that the number of wavelet coefficients we want to use is legal
#

if (is.na(IsPowerOfTwo(Ncoefs+1)))
	stop("Ncoefs has to be a power of two minus 1. E.g. =3, 7, 15, etc")

#
# Generate random starting coefficients, in parallel if possible
#

startpar <-  as.list(rep(2*Ncoefs, Nsims))	# Starting coefficients
startpar <- lapplyfn(startpar, irng)

#
# Apply optimizer to each set of starting coefficients using the
# objective fn lctsfn
#

ans <- lapplyfn(startpar, optim, fn=lctsfn, tsx=sf*tsx, tsy=sf*tsy,
	plot.it=plot.it, filter.number=filter.number, family=family,
	spec.filter.number=spec.filter.number, spec.family=spec.family,
	control=optim.control, ...)

#
# Get and store starting, ending coefficients, convergence info, minimum
# variance and pvalue for aech solution into a matrix or vector
#

startpar <- matrix(unlist(startpar), nrow=Nsims, byrow=TRUE)
endpar <- matrix(unlist(endparL <- lapplyfn(ans, getElement, name="par")),
	nrow=Nsims, byrow=TRUE)
convergence <- unlist(lapplyfn(ans, getElement, name="convergence"))
minvar <- unlist(lapplyfn(ans, getElement, name="value"))

pvals <- unlist(lapplyfn(endparL, getpvals, prodcomb.fn=prodcomb.fn,
	tsx=tsx, tsy=tsy, filter.number=filter.number, family=family,
	verbose=verbose, tos=BootTOS, Bsims=Bsims, lapplyfn=lapplyfn))

#
# Build return object and return it
#

l <- list(startpar=startpar, endpar=endpar, convergence=convergence,
	minvar=minvar, pvals=pvals, tsx=tsx, tsy=tsy,
	tsxname=deparse(substitute(tsx)),
	tsyname=deparse(substitute(tsy)),
	filter.number=filter.number, family=family,
	spec.filter.number=spec.filter.number, spec.family=spec.family)
class(l) <- "csFSS"
return(l)
}
