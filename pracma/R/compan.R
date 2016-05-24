###
### COMPAN.R  Polynom
###

compan <- function(p) {
	if (length(p) == 0) return(c())
	if ( !is.vector(p, mode="numeric"))
		stop("Argument p must be a vector of real numbers.")
	while(p[1] == 0 && length(p) >= 2) {
		p <- p[2:length(p)]
	}

	n <- length(p)
	if (n <= 1) {
		a <- c()
	} else {
		if (n == 2) {
			a <- -p[2]/p[1]
		} else {
			a <- diag(0, n-1, n-1)
			for (i in 2:(n-1)) {
				a[i, i-1] <- 1
			}
			a[1, ] <- -p[2:n]/p[1]
		}
	}
	return(a)
}
