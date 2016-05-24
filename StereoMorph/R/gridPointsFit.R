gridPointsFit <- function(p, nx, ny=NULL){	
	
	if(is.null(ny)) return(seq(0, length=nx, by=p[2]) + p[1])

	ny_add <- matrix(t(matrix(seq(0, length=ny, by=p[3]), nrow=ny, ncol=nx)), ncol=1)
	
	return(ny_add + rep(seq(0, length=nx, by=p[2]) + p[1], ny))
}
