find_coords <-
function(wgts,nns,N,n,m)
{
	W <- wgts
	M <- t(diag(1,N)-W)%*%(diag(1,N)-W)
	
	#calculate the eigenvalues and -vectors of M
	e <- eigen(M)
	
	#choose the eigenvectors belonging to the 
	#m smallest not-null eigenvalues
	#and re-scale data
	Y <- e$vectors[,c((N-m):(N-1))]*sqrt(N)
	
	return(Y)
}

