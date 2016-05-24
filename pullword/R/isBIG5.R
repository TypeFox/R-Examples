
##' Indicate whether the encoding of input string is BIG5.
##' 
##' @title Indicate whether the encoding of input string is BIG5.
##' @param string A character vector.
##' @param combine Whether to combine all the strings.
##' @return Logical value.
##' @author Jian Li <\email{rweibo@@sina.com}>

isBIG5 <- function(string, combine = FALSE)
{
	string <- .verifyChar(string)
	if (length(string)  == 1) {
		OUT <- .C("CWrapper_encoding_isbig5", 
				characters = as.character(string),  
				numres = 2L,PACKAGE='pullword')
		OUT <- as.logical(OUT$numres)
	} else {
		if (combine) {
			OUT <- isBIG5(paste(string, collapse = ""))
		} else {
			OUT <- as.vector(sapply(string, isBIG5))
		}
	}
	return(OUT)
}

