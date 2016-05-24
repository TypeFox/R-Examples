## ================================================
## Retrieve the color palette from a sequence object
## ================================================

cpal <- function(seqdata) {
	palette <- attr(seqdata,"cpal")
	return(palette)
}

