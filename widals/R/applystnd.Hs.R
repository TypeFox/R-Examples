applystnd.Hs <-
function( Hs0, x ) {

	sHs0 <- t( ( t(Hs0) - x$h.mean ) / x$h.sd )
	if( x$intercept ) { sHs0[ , 1] <- 1 / sqrt(x$n) }

	return(sHs0)
}
