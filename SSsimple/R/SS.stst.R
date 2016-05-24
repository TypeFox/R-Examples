SS.stst <-
function( F, H, Q, R, P0, epsilon, verbosity=0 ) {
	d <- ncol(H)
	n <- nrow(H)

	params <- .internal.chk.mod.params( F, H, Q, R, P0=P0, d=d, n=n )	
	F <- params$F ; H <- params$H ; Q <- params$Q ; R <- params$R ; P0 <- params$P0
	

	P.apos <- P0
	P.apri <- P0
	I <- diag(1, d)
	j <- 0
	stst <- FALSE
	while( !stst ) {

		j <- j + 1
		if( verbosity>0 ) { cat( paste("iteration: ", j) , "\n" ) }
		
		P.apri.temp <- F %*% P.apos %*% t(F) + Q
		##delta.apri <- abs(  sum( diag(P.apri) ) - sum( diag(P.apri.temp) )  )
        delta.apri <- sum( diag(P.apri-P.apri.temp)^2  )
		P.apri <- P.apri.temp
		
		K <- P.apri %*% t(H) %*% solve( H %*% P.apri %*% t(H) + R, tol=0 )

		P.apos.temp <- ( I - K %*% H ) %*% P.apri
        ##delta.apos <- abs(  sum( diag(P.apos) ) - sum( diag(P.apos.temp) )  )
		delta.apos <- sum( diag(P.apos - P.apos.temp)^2  )
		P.apos <- P.apos.temp
		
		if( delta.apri < epsilon   &  delta.apos < epsilon ) { stst <- TRUE }
			
        if (verbosity > 1) {
            print(P.apri) ; print(P.apos)
            print(delta.apri) ; print(delta.apos)
        }
		
	}
	cat( paste("coverged at iteration:", j) , "\n" )
	return( list( P.apri=P.apri, P.apos=P.apos ) )
}

