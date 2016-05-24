##
##  p a s c a l . R  Pascal Triangle
##


pascal <- function(n, k=0) {
    # Tests: P0 == P1 %*% t(P1);  P2 %^% 3 == eye
	stopifnot(is.numeric(n), length(n) == 1,
	          is.numeric(k), length(k) == 1)
	if (!(k %in% c(0, 1, 2)))
		stop("Argument 'k' must be 0, 1 or 2.")
	if (floor(n) != ceiling(n)) n <- floor(n)
	if (n <= 0) return(c())
	if (n == 1) return(c(1))

	p <- matrix(0, nrow=n, ncol=n)
	p[1, ] <- rep(1, n)
    for (i in 2:n) {
		p[i, ] <- cumsum(p[i-1, ])
	}
	if (k == 0) return(p)
	# k == 1:
	p1 <- matrix(0, nrow=n, ncol=n)
	for (j in 1:n) {
		p1[j, 1:j] <- diag(p[1:j, j:1]) * (-1)^(0:(j-1))
	}
	if (k == 2) {
		p1 <- rot90(p1, -1)
		if (n %% 2 == 0) p1 <- -p1
	}	
	return(p1)
}
