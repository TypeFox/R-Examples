#################################################
# helper function Wald test
WaldTest <- function( delta , vcov , R , nobs , cvec = NULL , eps=1E-10 ){
		delta.ii <- delta
		var.delta.ii <- vcov
		if ( is.null(cvec) ){
			cvec <- rep( 0 , nrow(R) )
					}
		Rdii <- R %*% delta.ii - cvec
		v1 <- R %*% var.delta.ii %*% t(R)
		diag(v1) <- diag(v1)+ eps	
		stat <- ( t( Rdii ) %*% solve( v1 ) %*% Rdii )[1,1]
		stats <- list()
		stats["X2"] <- stat
		stats["df"] <- nrow(R)
		stats["p"] <- 1 - stats::pchisq( stat , df = nrow(R) )
		l1 <- stats$X2 / stats$df - 1
		l1 <- ifelse( l1 < 0 , 0 , l1 / ( sum( nobs - 1 ) ) )
		stats["RMSEA"] <- sqrt(l1)
		return(stats)
				}
##############################################				