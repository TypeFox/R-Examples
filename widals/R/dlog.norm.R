dlog.norm <-
function(n, center, sd) {
	return( exp( rnorm(n, log( center ), sd) ) )
}
