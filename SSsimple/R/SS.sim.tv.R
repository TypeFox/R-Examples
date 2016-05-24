SS.sim.tv <-
function( F, H, Q, R, length.out, beta0=0 ) {
	if( class(H) == "list" ) { d <- ncol(H[[1]]) ; n <- nrow(H[[1]]) } else { d <- ncol(H) ; n <- nrow(H) }
	params <- .internal.chk.mod.params.short( P0=NULL, beta0, d, n )
	beta0 <- params$beta0

	T <- length.out
	Beta <- matrix(NA, T, d)
	Z <- matrix(NA, T, n)
	Y <- matrix(NA, T, n)

	for(j in 1:T) {
		if( class(F) == "list" ) { F.t <- F[[j]] } else { F.t <- F }
		if( class(H) == "list" ) { H.t <- H[[j]] } else { H.t <- H }
		if( class(Q) == "list" ) { Q.t <- Q[[j]] } else { Q.t <- Q }
		if( class(R) == "list" ) { R.t <- R[[j]] } else { R.t <- R }
########################### beta0 MUST BE VECTOR OF LENGTH D
		params <- .internal.chk.mod.params( F.t, H.t, Q.t, R.t, P0=NULL, beta0=NULL, d, n )	

		F.t <- params$F ; H.t <- params$H ; Q.t <- params$Q ; R.t <- params$R
		
		eta <- t( rmvnorm( 1, rep(0, d), Q.t ) )
		epsilon <- t( rmvnorm( 1, rep(0, n), R.t ) )
		if( j == 1 ) {
			Beta[ j, ] <- F.t %*% beta0 + eta
		} else {
			Beta[ j, ] <- F.t %*% Beta[ j-1, ] + eta
		}
		Y[ j, ] <- H.t %*% Beta[ j, ]
		Z[ j, ] <- Y[ j, ] + epsilon
	}
	return( list( Beta=Beta, Y=Y, Z=Z ) )
}

