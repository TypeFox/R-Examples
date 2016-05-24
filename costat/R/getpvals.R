getpvals <-
function(par, prodcomb.fn, tsx, tsy, filter.number,
	family=c("DaubExPhase", "DaubLeAsymm"),
	verbose=FALSE, tos=BootTOS, Bsims=100, lapplyfn=lapply)
{
family <- match.arg(family)
#
# Compute p-value for a particular linear combination
#

#
# Build the combination on the the two time series: Z_t
#

ans1 <- prodcomb.fn(par,tsx=tsx,tsy=tsy, filter.number=filter.number,
	family=family)

#
# Compute the appropriate test of stationarity
#
ans2 <- tos(ans1, Bsims=Bsims, lapplyfn=lapplyfn)
#
# Compute the p-value of the test
#
pval <- plotBS(ans2$Bootvals, plot=FALSE, verbose=verbose)
#
# Return it
#
return(pval)
}
