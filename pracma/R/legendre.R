##
##  l e g e n d r e . R  Legendre Functions
##


legendre <- function(n, x) {
	stopifnot(is.numeric(x), is.numeric(n),
	          length(n) == 1, floor(n) == ceiling(n), n >= 0)
	x <- c(x)
	N <- length(x)

	if (n == 0) return(rep(1, 10))

    # generate the Legendre polynomials up to degree n
    Lp <- matrix(0, n+1, n+1)
    Lp[1, n+1] <- Lp[2, n] <- 1
    if (n > 1) {
        for (i in 3:(n+1)) {
	        j <- i-1
	        Lp[i, (n-i+2):(n+1)] <- (2*j-1)/j * c(Lp[i-1, (n-i+3):(n+1)], 0) -
	                      (j-1)/j * Lp[i-2, (n-i+2):(n+1)]
        }
    }
    lp <- Lp[n+1, ]

    # associated Legendre functions up to order n
    L <- matrix(NA, n+1, N)
    L[1, ] <- polyval(lp, x)

    for (j in 1:n) {
	    lp <- polyder(lp)
	    L[j+1, ] <- (-1)^j * sqrt((1-x^2)^j) * polyval(lp, x)
    }

    return(L)
}
