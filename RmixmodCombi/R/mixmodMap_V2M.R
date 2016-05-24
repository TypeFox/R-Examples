mixmodMap_V2M <-
function(z, n = length(z), K = max(z))
{
	A = z + cumsum(c(0,rep(K,(n-1))))
	B = numeric(n * K)
	B[A] = 1
	M <- matrix(B, nrow = n, ncol = K, byrow = T)
	M
}
