mixmodMap <-
function(tau, n = nrow(tau), K = ncol(tau))
{
	M <- matrix(0, nrow = n, ncol = K)
	tauMax <- apply(tau, 1, which.max)
	A = tauMax + cumsum(c(0,rep(K,(n-1))))
	B = numeric(n * K)
	B[A] = 1
	M <- matrix(B, nrow = n, ncol = K, byrow = T)
	M
}
