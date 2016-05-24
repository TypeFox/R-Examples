func <- function(v,fu,m,nfmax,da,ga,tr,nql)
{
	for(i in 1:m){
		ii <- i*(i+1)/2
		if(v[ii] != 0){
			fu[ii] <- v[ii]+da*((1-ga)*abs(v[ii])+ga*tr)
		}else{
			fu[ii] <- da*ga*tr
		}
	}

	dchole <- .Fortran("dchole",fu=as.double(fu),as.integer(m),as.integer(nql),idpos=as.integer(0),PACKAGE="marqLevAlg")
	idpos <- dchole$idpos
	fu <-dchole$fu
	result <- list(idpos=idpos,fu=fu)
	return(result)
}