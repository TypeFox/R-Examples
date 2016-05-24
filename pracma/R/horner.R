##
##  h o r n e r . R  Horner Scheme
##


# Horner's rule to compute p(x) and p'(x) vectorized for the
# polynomial  p = p_1*x^n + p_2*x^{n-1} + ... + p_n*x + p_{n+1}
horner <- function(p, x) {
    if (length(p) == 0 || length(x) == 0) return(NULL)

	n <- length(p); m <- length(x)
	if (n == 0) { y <- dy <- rep(NA, m) }
	else if (n == 1) { y <- rep(p, m); dy <- rep(0, m) }
	else {
		y <- p[1]; dy <- 0
		for (i in 2:n) {
			dy <- dy * x + y
			y  <-  y * x + p[i]
		}
	}
	return(list(y = y, dy = dy))
}


# Deflated Horner scheme that returns p(x) and the polynomial q with
# p(x) = q(x) * (x - x0) + r, r constant, and r = 0 if x0 is a root of p.
hornerdefl <- function(p, x) {
    if (length(p) == 0 || length(x) == 0) return(NULL)

    n <- length(p) -1                   # degree of polynomial
    q <- numeric(n+1)
    q[1] <- p[1]
    for (j in 2:(n+1))
        q[j] <- p[j] + q[j-1]*x
    return(list(y = q[n+1], q = q[1:n]))
}
