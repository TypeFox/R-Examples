
.nchar <- function(x) {
	x[is.na(x)] <- ''
	sapply(strsplit(x, NULL), length)
}


