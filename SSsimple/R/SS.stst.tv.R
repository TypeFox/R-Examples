SS.stst.tv <-
function( F, H, Q, R, P0, epsilon, verbosity=0 ) {
	if( class(H) == "list" ) { d <- ncol(H[[1]]) ; n <- nrow(H[[1]]) } else { d <- ncol(H) ; n <- nrow(H) }
	params <- .internal.chk.mod.params.short( P0, beta0=NULL, d, n )
	P0 <- params$P0
	
	
	P.apos <- P0
	P.apri <- P0
	I <- diag(1, d)
	j <- 0
	stst <- FALSE

	while( !stst ) {

		j <- j + 1
		
		if( verbosity>0 ) { cat( paste("iteration: ", j) , "\n" ) }
		if( class(F) == "list" ) { F.t <- F[[j]] } else { F.t <- F }
		if( class(H) == "list" ) { H.t <- H[[j]] } else { H.t <- H }
		if( class(Q) == "list" ) { Q.t <- Q[[j]] } else { Q.t <- Q }
		if( class(R) == "list" ) { R.t <- R[[j]] } else { R.t <- R }
		params <- .internal.chk.mod.params( F.t, H.t, Q.t, R.t, P0=NULL, beta0=NULL, d, n )
		F.t <- params$F ; H.t <- params$H ; Q.t <- params$Q ; R.t <- params$R

		
		P.apri.temp <- F.t %*% P.apos %*% t(F.t) + Q.t
		delta.apri <- sum( diag(P.apri - P.apri.temp)^2  )
		P.apri <- P.apri.temp
		
		K <- P.apri %*% t(H.t) %*% solve( H.t %*% P.apri %*% t(H.t) + R.t, tol=0 )
		
		P.apos.temp <- ( I - K %*% H.t ) %*% P.apri
		delta.apos <- sum( diag(P.apos - P.apos.temp)^2  )
		P.apos <- P.apos.temp
		
		if( delta.apri < epsilon   &  delta.apos < epsilon ) { stst <- TRUE }
		
		if( verbosity>1 ) { print( P.apri ) ; print( P.apos ) }
		if( verbosity>1 ) { print( delta.apri ) ; print( delta.apos ) }
		
	}
	cat( paste("coverged at iteration:", j) , "\n" )
	return( list( P.apri=P.apri, P.apos=P.apos ) )
}

