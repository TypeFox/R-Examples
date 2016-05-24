"reglin" <-
function (x, y=NULL, xmin=min(x), n=length(x), deltat=(max(x)-min(x))/(n-1), rule=1) {
	# We use approx() for the calculations
	# but we first need to construct xout
	if (n <= 0) 
    	stop("reglin requires n >= 1")
    if (deltat <= 0) 
    	stop("reglin requires deltat > 0")
	xout <- 0:(n-1) * deltat + xmin
	# Make sure data are sorted in increasing order according to x
	srt <- sort.list(x)
	x <- x[srt]
	y <- y[srt]
	res <- approx(x, y, xout, method="linear", rule=rule)
	res
}
