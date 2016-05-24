stnd.Hst.ls <-
function( Hst.ls, Hst0.ls=NULL ) {
    tau <- length(Hst.ls)
	big.sum <- 0
	for(i in 1:tau) { 
		big.sum <- big.sum +  apply( Hst.ls[[i]] , 2, mean )
	}
	h.mean <- big.sum / tau
	sHst.ls <- list()
	big.sum.mx <- 0
	for(i in 1:tau) { 
		sHst.ls[[i]] <- t( t(Hst.ls[[i]]) - h.mean )
		big.sum.mx <- big.sum.mx + crossprod( sHst.ls[[i]] )
	}
	cov.mx <- big.sum.mx / tau
	sqrtXX <- 1 / sqrt( diag(cov.mx) )
	for(i in 1:tau) { 
		sHst.ls[[i]] <- t( t(sHst.ls[[i]]) * sqrtXX )
	}
	
	sHst0.ls <- NULL
	if( !is.null(Hst0.ls) ) {
		sHst0.ls <- list()
		for(i in 1:tau) { 
			sHst0.ls[[i]] <- t( ( t(Hst0.ls[[i]]) - h.mean ) * sqrtXX )
		}
	}
	ls.out <- list( "sHst.ls"=sHst.ls, "sHst0.ls"=sHst0.ls, "h.mean"=h.mean, "h.sd"= 1/sqrtXX )
    return( ls.out )
}
