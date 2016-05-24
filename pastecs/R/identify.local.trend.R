"identify.local.trend" <-
function(x, ...) {
	x2 <- x
	class(x2) <- NULL
	x2 <- as.ts(x2)
	pos <- identify(x2, ...)
	# Calculate trends for the successive segments using slope and k
	n <- length(pos)
	k <- attr(x, "k")
	if (n > 1) {
		x2.val <- x2[pos]
		slopes <- (x2.val[2:n] - x2.val[1:(n-1)]) / (pos[2:n] - pos[1:(n-1)]) + k
	} else {
		slopes <- NA
	}
	res <- list(pos=pos, trends=slopes, k=k)
	res
}
