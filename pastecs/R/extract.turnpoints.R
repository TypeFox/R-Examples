"extract.turnpoints" <-
function(e, n, no.tp=0, peak=1, pit=-1, ...) {
	if (missing(n)) n <- -1
	res <- rep(no.tp, length.out=e$n)
	res[e$pos[e$peaks]] <- peak
	res[e$pos[e$pits]] <- pit
	# Keep only the first n points
	if (n < length(res) & n > 0) res <- res[1:n]
	res
}
