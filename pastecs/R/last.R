"last" <-
function(x, na.rm=FALSE) {
	if (na.rm) 
		x <- x[!is.na(x)]
	x[length(x)]
}
