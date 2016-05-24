##
##  p a d e . R  Pade Approximation
##


pade <- function(p1, p2 = c(1), d1 = 5, d2 = 5) {
	stopifnot(is.numeric(p1), is.numeric(p2),
	    is.numeric(d1), length(d1) == 1, floor(d1) == ceiling(d1), d1 >= 0,
	    is.numeric(d2), length(d2) == 1, floor(d2) == ceiling(d2), d2 >= 0)

	if (d1 == 0 && d2 == 0)
	    return(list(r1 = p1[length(p1)], r2 = p2[length(p2)]))

	z   <- rep(0, d1 + d2 + 3)
	p2 <- rev(c(z, p2))[1:(d1+d2+3)]
	p1 <- rev(c(z, p1))[1:(d1+d2+3)]

	L <- Toeplitz(p2[1:(d1+d2)]); L[upper.tri(L)] <- 0
	R <- Toeplitz(p1[1:(d1+d2)]); R[upper.tri(R)] <- 0

    # generate the linear system of coefficient equations
    D1 <- if (d1 > 0) 1:d1 else c()
    D2 <- if (d2 > 0) 1:d2 else c()
	A <- cbind(L[, D1], -R[, D2])
	b <- p2[1]*p1[2:(d1+d2+1)] - p1[1]*p2[2:(d1+d2+1)]

	# imitate pinv() if some eigenvalues are zero
	P <- svd(A); U <- P$u; V <- P$v; s <- P$d
	r <- sum(s > max(dim(A)) * max(s) * .Machine$double.eps)
	if (r == 0) {
	    pinvA <- matrix(0, nrow = nrow(A), ncol = ncol(A))
	} else {
	    S <- diag(1/s[1:r])
	    pinvA <- V[, 1:r] %*% S %*% t(U[, 1:r])
    }

	# solve the linear system of coefficient equations
	B <- zapsmall(pinvA %*% b)

	# reconstruct the rational function
	r1 <- rev(c(p1[1], B[1:d1]))
	r2 <- rev(c(p2[1], B[(d1+1):length(B)]))

	# scale such that max(r2) = 1
	rmax <- max(abs(r2))
	if (rmax == 0) rmax <- 1
	r1 <- r1 / rmax
	r2 <- r2 / rmax

	# return numerator and denominator
	return(list(r1 = r1, r2 = r2))
}
