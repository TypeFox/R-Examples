fixdmat <- function(v)
{
# Converts a row-wise distance vector to a full distance matrix.
	n <- (1 + sqrt(1 + 8 * length(v)))/2
	dist.m <- matrix(0, n, n)
	dist.m[row(dist.m) < col(dist.m)] <- v
	dist.m <- dist.m + t(dist.m)
	dist.m
}
