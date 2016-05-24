amari <-
function(V, W, orth = FALSE)
{
### a metric between two orthonormal matrices
	# V and W should be square
	if(orth) A <- abs(t(V) %*% W) else A <- abs(solve(V, W))
	rmax <- apply(A, 1, max)
	cmax <- apply(A, 2, max)
	rsum <- apply(A, 1, sum)
	csum <- apply(A, 2, sum)
	(sum(rsum/rmax - 1) + sum(csum/cmax - 1))/(2 * nrow(A))
}

