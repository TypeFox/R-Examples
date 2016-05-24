# sequential Mann-Kendall test on rank time series (after Sneyers 1990)
# detects approximate potential trend turning points in time series
# returns the progressive and retrograde time series of Kendall's
# normalized tau.
# input: vector x of the time series, length n
# output : progressive/retrograde series, length n-1 + NA at the beginning/end,
# and a vector of indices of the original times where potential approximate
# trend turning points are situated
seqMK <- function(x) {
	n <- length(x)
	y <- rev(x)
	prog <- retr <- vector("numeric",n)
	tp <- vector("logical",n)
	prog[1] <- retr[1] <- NA
	tp[1] <- tp[n] <- FALSE
	if(n < 2)
	    stop("seqMK: not enough finite observations")
	# progressive and retrograde series
	for (i in 2:n) {
		prog[i] <- tau(x[1:i])
		retr[i] <- tau(y[1:i])
	}
	retr <- rev(retr)
	diff <- prog-retr
	# index vector of crossing points
	for (i in 2:(n-2)) {
		if(sign(diff[i])==sign(diff[i+1])) { tp[i+1] <- FALSE }
		else { tp[i+1] <- TRUE }
	}
	return(list(prog=prog, retr=retr, tp=tp))
}
