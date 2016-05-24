stnd.Ht <-
function(Ht, n) {
	h.mean <- apply(Ht, 2, mean)
	sHt <- t( t(Ht) - h.mean )
	sHt <- t( t(sHt) / apply( sHt, 2, function(x) { sqrt( sum(x^2) ) } ) )
	sHt <- sHt * sqrt( nrow(Ht) / n )
	return(sHt)
}
