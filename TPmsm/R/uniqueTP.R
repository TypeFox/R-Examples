uniqueTIME <- function(object, s, t) {
	return( .Call(Rf_uniqueTIME, object, s, t, PACKAGE="TPmsm") )
}

uniqueCOV <- function(object, x) {
	return( .Call(Rf_uniqueCOV, object, x, PACKAGE="TPmsm") )
}
