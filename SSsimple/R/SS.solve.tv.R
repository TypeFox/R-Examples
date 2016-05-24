SS.solve.tv <-
function( Z, F, H, Q, R, length.out, P0, beta0 ) {
	if( class(H) == "list" ) { d <- ncol(H[[1]]) ; n <- nrow(H[[1]]) } else { d <- ncol(H) ; n <- nrow(H) }
	params <- .internal.chk.mod.params.short( P0, beta0, d, n )
	P0 <- params$P0 ; beta0 <- params$beta0

	T <- length.out
	
	B.apri <- matrix(NA, T, d)
	B.apos <- matrix(NA, T, d)
	Z.hat.apri <- matrix(NA, T, n)
	Z.hat.apos <- matrix(NA, T, n)
	P <- P0
	I <- diag(1, d)
	
	for(j in 1:T) {
		if( class(F) == "list" ) { F.t <- F[[j]] } else { F.t <- F }
		if( class(H) == "list" ) { H.t <- H[[j]] } else { H.t <- H }
		if( class(Q) == "list" ) { Q.t <- Q[[j]] } else { Q.t <- Q }
		if( class(R) == "list" ) { R.t <- R[[j]] } else { R.t <- R }
		params <- .internal.chk.mod.params( F.t, H.t, Q.t, R.t, P0=NULL, beta0=NULL, d, n )
		F.t <- params$F ; H.t <- params$H ; Q.t <- params$Q ; R.t <- params$R

		if(j == 1) {
			B.apri[ j, ] = F.t %*% beta0
			
		} else {
			B.apri[ j, ] = F.t %*% B.apos[ j-1, ]
		}
		Z.hat.apri[j, ] <- H.t %*% B.apri[ j, ]
	
		P <- F.t %*% tcrossprod( P, F.t ) + Q.t
		
		PtH <- tcrossprod( P, H.t )
		
		K <- PtH %*% solve( H.t %*% PtH + R.t, tol=0 )
		
		B.apos[ j, ] <- B.apri[ j, ] + K %*% ( Z[j, ] - H.t %*% B.apri[ j, ] )
		Z.hat.apos[j, ] <- H.t %*% B.apos[ j, ]
		
		P <- ( I - K %*% H.t ) %*% P
		
		
	}
	
	return( list( B.apri=B.apri, B.apos=B.apos, fitted.apri=Z.hat.apri, fitted.apos=Z.hat.apos ) )
}

