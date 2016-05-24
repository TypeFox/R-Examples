stnd.Hs <-
function(Hs, Hs0=NULL, intercept=TRUE) {
	n <- nrow(Hs)
	h.mean <- apply(Hs, 2, mean)
	h.sd <- apply( t(t(Hs)-h.mean) , 2, function(x) { sqrt( sum(x^2) ) } )
	h.sd[ h.sd == 0 ] <- 1
	sHs <- t( ( t(Hs) - h.mean ) / h.sd )
	if( intercept ) { sHs[ , 1] <- 1 / sqrt(n) }
	
	sHs0 <- NULL
	if( !is.null(Hs0) ) {
		sHs0 <- t( ( t(Hs0) - h.mean ) / h.sd )
		if( intercept ) { sHs0[ , 1] <- 1 / sqrt(n) }
	}
	
	ls.out <- list( "sHs"=sHs, "sHs0"=sHs0, "h.mean"=h.mean, "h.sd"=h.sd, "n"=n, "intercept"=intercept )
	return(ls.out)
}
