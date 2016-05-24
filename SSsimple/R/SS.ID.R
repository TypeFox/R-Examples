SS.ID <-
function( Z, d, rsN = NULL ) {
	T <- nrow(Z)
	n <- ncol(Z)
	time.t <- as.integer( floor(T/2) )
	if( is.null(rsN) ) { 
		r <- min(10, as.integer( floor( time.t / n ) - 3 ) )
		s <- min(20, as.integer( floor( time.t / n ) - 2 ) ) 
		N <- as.integer( floor(T/2) - max(r,s) - 2 )
	} else {
		r <- rsN[1] ; s <- rsN[2] ; N <- rsN[3]
	}
	
	YY <- matrix( NA, n*r,N )
	for( tt in time.t:(time.t+N-1) ) { YY[ , tt - time.t + 1 ]  <- as.vector( t( Z[ tt:(tt+r-1), ] ) ) }
	
	PP <- matrix( NA, n*s, N )
	for( tt in time.t:(time.t+N-1) ) { PP[ , tt - time.t + 1 ] <- as.vector( t( Z[ (tt-1):(tt-s), ] ) ) }
	
	GG <- tcrossprod( YY , PP ) / N
	
	W1 <- diag(1, r*n)
	W2 <- solve( tcrossprod(PP) / N ) %*% PP ##### N used in next step.  Not included here so we can use W2 later
	
	GG.hat <- W1 %*% GG %*% W2
	
	GGhs <- svd( GG.hat )
	
	U1 <- GGhs$u[ , 1:d, drop=FALSE ]
	RR <- diag(1, d)
	
	O.hat <- solve( W1 ) %*% U1 %*% RR
	
	H.hat <- O.hat[1:n , 1:d, drop=FALSE]
	
	O1 <- O.hat[ (n+1):(n*r) , 1:d, drop=FALSE ]
	O2 <- O.hat[ (1):(n*(r-1)) , 1:d, drop=FALSE ]
	
	F.hat <- tcrossprod( solve( crossprod(O2) ) , O2 ) %*% O1
	
	LL <- tcrossprod( solve( RR ) , U1 )
	
	YY.hat <- tcrossprod( YY , PP ) %*% solve( tcrossprod(PP) ) %*% PP
	
	XX.hat <- LL %*% YY.hat
	
	nu.mx <- matrix(NA, d, ncol(XX.hat)-1)
	for(i in 2:ncol(XX.hat)) {
		nu.mx[ , i-1] <- XX.hat[ , i ] - F.hat %*% XX.hat[ , i-1 ]
	}
	Q.hat <- var( t(nu.mx) )
	
	ve.mx <- matrix(NA, n, ncol(XX.hat))
	for(i in 1:ncol(XX.hat)) {
		ve.mx[ , i] <- Z[ time.t+i-2, ] - H.hat %*% XX.hat[ , i ] ########### OK
	}
	R.hat <- var( t(ve.mx) )
	
	return( list( F.hat=F.hat, H.hat=H.hat, Q.hat=Q.hat, R.hat=R.hat ) )
	
}

