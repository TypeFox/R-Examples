
.verifyChar <- function(vec) {
	if (!is.atomic(vec)) stop("Must be an atomic vector!", call. = FALSE)
	OUT <- as.character(vec)
	return(OUT)
}

.verifyNum <- function(vec) {
	if (!is.atomic(vec)) stop("Must be an atomic vector!")
	OUT <- suppressWarnings(as.numeric(vec))
	if (length(OUT[!is.na(OUT)]) == 0) stop("Does not contain any numeric elements!", call. = FALSE)
	return(OUT)
}

