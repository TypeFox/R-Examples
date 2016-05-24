##
##  s q r t m . R  Matrix Square and p-th Roots
##


sqrtm <- function(A, kmax = 20, tol = .Machine$double.eps^(1/2)) {
	stopifnot(is.numeric(A), is.matrix(A))
	if (nrow(A) != ncol(A))
	    stop("Matrix 'A' must be square.")

	# should be "try(solve(A))"
	P0 <- A; Q0 <- diag(nrow(A))

	k <- 0  # then k <- 1
	while (norm(A - P0 %*% P0, 'F') > tol && k < kmax) {
		P1 <- 0.5 * (P0 + solve(Q0))
		Q1 <- 0.5 * (Q0 + solve(P0))
		P0 <- P1
		Q0 <- Q1
		k <- k + 1
	}

	# k < 0 if iteration has not converged.
	if (k >= kmax) k <- -1

	# return sqrtm(A) and sqrtm(A)^-1
	return(list(B = P0, Binv = Q0, k = k, acc = norm(A - P0 %*% P0, 'F')))
}


signm <- function(A, kmax = 20, tol = .Machine$double.eps^(1/2)) {
    A %*% sqrtm(A %*% A)$Binv
}


rootm <- function(A, p, kmax = 20, tol = .Machine$double.eps^(1/2)) {
	stopifnot(is.numeric(A), is.matrix(A))
	if (nrow(A) != ncol(A))
	    stop("Matrix 'A' must be square.")

	n <- nrow(A)
	A0 <- A
	p0 <- p

    # err <- try(solve(A), silent = TRUE)
    # if (class(err == "try-error")) ... else ...

	if (p %% 2 == 1) {
		A <- A %*% A
		p <- 2*p
	} else {
		while (p %% 4 == 0) {
			A <- sqrtm(A)$B
			p <- p / 2
		}
	}

	N <- N0 <- 1
	acc <- Inf
	k <- 0
	while (acc > tol && k < kmax) {  # && acc <= accp
		N <- 2*N
		wN <- cos(2*pi/N) + 1i * sin(2*pi/N)  # N-th root of unity

		# summing with roots of unity
		S <- solve(A)/4
		for (j in 1:(N-1)) {
			B <- solve( A - ((1-wN^j)/(1+wN^j))^p * diag(n) )
			S <- S + B * wN^j/(1+wN^j)^2
		}
		S <- 2*p*sin(pi/p)/N * A %*% S
		S <- Re(S)

		Sp <- S
		for (j in 1:(p0-1)) Sp <- Sp %*% S

		acc <- norm(A0 - Sp, 'F')
		k <- k + 1
	}

	# k < 0 if iteration has not converged.
	if (k >= kmax) k <- -1

	return(list(B = S, k = k, acc = acc))
}


