#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 17:15:39 +0100 (ven, 02 dic 2005) $
lyap_k <- function(series, m, d, t, k=1, ref, s, eps) {
	checkEmbParms(series, m, d, t, s, ref)
	series <- as.ts(series)
	n <- length(series) - (m-1)*d - s 
	if(ref<0) ref <- n
	trash <- numeric()
	ref <- 1:ref
	cat("Finding nearests\n")
	nearest <- find_knearests(series, m=m, d=d, t=t, ref=length(ref), s=s, eps=eps, k=k)
	trash <- apply(nearest, 1, function(x) any(is.na(x)))
	ref <- ref[!trash]
	if(length(ref)==0) stop("not enough neighbours found")
	cat("Keeping ", length(ref)," reference points\n")
	cat("Following points\n")
	res <- follow_points(series, m=m, d=d, s=s, ref=ref, nearest=nearest, k=k)
	ts(res, frequency = frequency(series), start=0)
}

lyap <- function(dsts, start, end) {
	dsts <- as.ts(dsts)
	xtsp <- tsp(dsts)
	seq(xtsp[1], xtsp[2], by = 1/xtsp[3])
	sf <- window(dsts, start, end)
	lambda <- tsp(sf)
	lambda <- seq(lambda[1], lambda[2], by=1/lambda[3])
	lm(sf~lambda, data=data.frame(sf = sf, lambda = lambda))$coeff
}
