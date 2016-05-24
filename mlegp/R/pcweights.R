'pcweights' <- function(Y, weights.num = NULL, cutoff = 99) {
	if (is.null(weights.num)) {
		weights.num = numSingularValues(Y, cutoff = cutoff)
	}

	s = svd(Y)
	ans = NULL
	ans$UD = s$u %*% diag(s$d)[,1:weights.num]
	ans$Vprime = t(s$v)[1:weights.num,]
	return (ans)
}
