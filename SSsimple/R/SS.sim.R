SS.sim <-
function( F, H, Q, R, length.out, beta0=0 ) {
	
	d <- ncol(H)
	n <- nrow(H)
	T <- length.out
	params <- .internal.chk.mod.params( F, H, Q, R, P0=NULL, beta0=beta0, d, n )	
	F <- params$F ; H <- params$H ; Q <- params$Q ; R <- params$R ; beta0 <- params$beta0
	Beta <- matrix(NA, T, d)
	Z <- matrix(NA, T, n)
	Y <- matrix(NA, T, n)

	Eta <- rmvnorm( T, rep(0, d), (Q) )

	Epsilon <- rmvnorm( T, rep(0, n), R )

	for(j in 1:T) {
		if( j == 1 ) {
			Beta[ j, ] <- F %*% beta0 + Eta[ j, ]
		} else {
			Beta[ j, ] <- F %*% Beta[ j-1, ] + Eta[ j, ]
		}
		Y[ j, ] <- H %*% Beta[ j, ]
		Z[ j, ] <- Y[ j, ] + Epsilon[ j, ]
	}
	return( list( Beta=Beta, Y=Y, Z=Z ) )
}

