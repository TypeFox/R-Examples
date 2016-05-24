full <- function(v)
{
# converts a lower-triangular distance vector to a symmetric matrix
	n <- (1 + sqrt(1 + 8 * length(v)))/2
        if(abs(n - round(n)) > 0.0000001)
		stop("Matrix not square.")
	n <- round(n)
        full <- matrix(0, n, n)
        full[lower.tri(full)] <- v
	full2 <- t(full)
	diag(full2) <- 0
        full + full2
}
