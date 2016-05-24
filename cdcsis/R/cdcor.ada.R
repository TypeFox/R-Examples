cdcor.ada <-
function(x, y, z, tol=1e-1,index=1){
	if(  is.null(dim(z)) == TRUE ) {
		width <- bw(x, y ,z,index)
		if (width<tol) width <- hpi(z)
	out <- cdcor(x, y, z, width,index)$mcdcor 
}

if(  is.null(dim(z)) == FALSE ) {
	width <- rep(NA, dim(z)[2])
	for( i in 1:dim(z)[2]) {
		width[i] <- bw(x, y, z[,i],index)
		if (width[i]<tol) width[i] <- hpi(z[,i])
	}
	w <- diag(width)
	out <- cdcor(x, y, z, w,index)$mcdcor 
}
return(list(cdcor=out,width=width))
}
