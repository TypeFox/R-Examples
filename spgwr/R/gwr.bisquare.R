gwr.bisquare <- function(dist2, d) {
	d2 <- d^2
	w <- ifelse(dist2 > d2, 0, (1 - (dist2/d2))^2)
	w
}
