kernel_twowayx <-
function(Z, n, p)
{
	K = matrix(0, n, n)	
	aux <- .C("kernel_twowayx", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), as.double(as.vector(K)))[[4]]
	matrix(aux, nrow=n)
}

