SS.solve <-
function( Z, F, H, Q, R, length.out, P0, beta0=0 ) {
	d <- ncol(H)
	n <- nrow(H)
	T <- length.out
	params <- .internal.chk.mod.params( F, H, Q, R, P0=P0, beta0=beta0, d, n )	
	F <- params$F ; H <- params$H ; Q <- params$Q ; R <- params$R ; P0 <- params$P0 ; beta0 <- params$beta0
	
	B.apri <- matrix(NA, T, d)

	B.apos <- matrix(NA, T, d)

	P <- P0
	I <- diag(1, d)

	for(j in 1:T) {
		
		if(j == 1) {
			B.apri[ j, ] = F %*% beta0

		} else {
			B.apri[ j, ] = F %*% B.apos[ j-1, ]
		}
		
		
		P <- F %*% tcrossprod( P, F ) + Q
		
		PtH <- tcrossprod( P, H )
		
		K <- PtH %*% solve( H %*% PtH + R, tol=0 )
		
		B.apos[ j, ] <- B.apri[ j, ] + K %*% ( Z[j, ] - H %*% B.apri[ j, ] )

		P <- ( I - K %*% H ) %*% P
		
		
	}
	
	return( list( B.apri=B.apri, B.apos=B.apos ) )
}

