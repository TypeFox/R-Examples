kernel_wIBS <-
function(Z, n, p, weights)
{
	#given_weight = 1;
	#if( is.null(weights)){
    #		weights = rep(0, p);
	#	given_weight = 0;
	#}
	# K = matrix(0, n, n)
	#aux = .C("kernel_wIBS", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), as.integer(given_weight),
	#as.double(weights), as.double(as.vector(K)))[[6]]

	K = diag(1, n, n)
	aux = .C("kernel_wIBS", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), 
		    as.double(weights), as.double(as.vector(K)))[[5]]
	matrix(aux, nrow=n)
}

