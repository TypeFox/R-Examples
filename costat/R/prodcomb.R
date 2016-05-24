prodcomb <-
function (cfs, tsx, tsy, filter.number = 1,
    family = c("DaubExPhase", "DaubLeAsymm"), 
    all = FALSE) 
{
family <- match.arg(family)
#
# Produce time-varying linear combination of tsx and tsy using coefficients
# in cfs
#
#
# Get dimensions of coefficients and time series
#
lcfs <- length(cfs)
lts <- length(tsx)
#
# Breakout alpha and beta coefficients
alpha <- cfs[1:(lcfs/2)]
beta <- cfs[(lcfs/2 + 1):lcfs]
#
# Turn wavelet coefficients into function
#
v <- coeftofn(alpha = alpha, beta = beta, n = lts,
	filter.number = filter.number, family = family)
#
# Form time-varying linear combination, Z_t
lcts <- v$alpha * tsx + v$beta * tsy
#
# Decide whether to return just the combination or everything.
#
if (all == FALSE) 
	return(lcts)
else {
	l <- list(lcts = lcts, alpha = v$alpha, beta = v$beta)
	return(l)
	}
}
