unif.mh <-
function(n, center, sd) {
	w <- sd * sqrt(3)
	a <- center - w
	b <- center + w
	x <- runif(n, a, b)
	return(x)
}
