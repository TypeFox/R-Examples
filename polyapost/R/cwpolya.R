
cwpolya <- function(samp, wts, k) {
	if (missing(wts)) wts <- rep(1, length(samp))
	if (length(samp) != length(wts))
		stop("length(samp) != length(wts)")
	if (any(wts < 0) || all(wts <= 0))
		stop("wts are foobar")
	k <- as.integer(k)
	if (k <= 0)
		stop("k must be positive integer")
       out<-.C(C_cwpolya,
       x=as.double(c(samp, rep(0, k))),
       w=as.double(cumsum(wts)),
       n=as.integer(length(samp)),
       k=as.integer(k))
       return(out$x)
}
