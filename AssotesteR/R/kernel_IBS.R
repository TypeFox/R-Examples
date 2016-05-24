kernel_IBS <-
function(Z, n, p)
{
	K = diag(1, n, n)	
	aux = .C("kernel_IBS", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), as.double(as.vector(K)))[[4]]
	matrix(aux, nrow=n)
}

