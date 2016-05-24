###-----------------------------------------------------------------------------
###
### getbynames - get list elements by name and if no such element exists, NA 
###
###

getbynames <- function (x, e) {
	x <- x [e]
	if (length (x) > 0) {
		if (is.character (e)) 
			names (x) <- e
		x [sapply (x, is.null)] <- NA
		x
	} else {
		list ()
	}
}
