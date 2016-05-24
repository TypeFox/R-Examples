##
##  r o r t h o . R  Generate Random Orthogonal Matrix
##


rortho <- function(n) {
	stopifnot(is.numeric(n), length(n) == 1, floor(n) == ceiling(n), n >= 1)
	if (n == 1) return(matrix(1, 1, 1))

    H <- diag(1, n, n)
    H[2:n, 2:n] <- H[2:n, 1 + sample.int(n-1,n-1)]

	u0 <- runif(n)
	u1 <- u0 + sqrt(sum(u0^2)) * c(1, rep(0, n-1))
	u1 <- u1/sqrt(sum(u1^2))

	Q <- diag(1, n, n) - 2 * as.matrix(u1) %*% t(u1)
	return(Q)
}
