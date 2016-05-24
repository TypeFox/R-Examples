TCS <-
function(basis, cov, lambda){
	cov[abs(cov)<lambda]=0
	smooth=basis%*%cov%*%t(basis)
	return(smooth)	
}
