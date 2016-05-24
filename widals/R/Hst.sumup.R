Hst.sumup <-
function(Hst.ls, Hs=NULL, Ht=NULL) {
	tau <- length(Hst.ls)
    if(tau < 1) { tau <- nrow(Ht) }
    if(is.null(tau)) { tau <- 10 ; cat("tau assumed to be 10.", "\n") }
    n <- nrow(Hst.ls[[1]])
    if(is.null(n)) { n <- nrow(Hs) }
	big.sum <- 0
	for(i in 1:tau) {
		if( !is.null(Ht) ) { Ht.mx <- matrix( Ht[i, ], n, ncol(Ht), byrow=TRUE ) } else { Ht.mx <- NULL }
		big.sum <- big.sum + crossprod( cbind( Hs, Ht.mx, Hst.ls[[i]] ) )
	}
	return( big.sum )
}
