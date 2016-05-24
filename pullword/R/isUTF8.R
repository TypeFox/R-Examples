
##' Indicate whether the encoding of input string is UTF-8.
##' 
##' @title Indicate whether the encoding of input string is UTF-8.
##' @param string A character vector.
##' @param combine Whether to combine all the strings.
##' @return Logical value.
##' @author Jian Li <\email{rweibo@@sina.com}>

isUTF8 <- function(string, combine = FALSE)
{
	string <- .verifyChar(string)
	if (length(string)  == 1) {
		OUT <- .C("CWrapper_encoding_isutf8", 
				characters = as.character(string),  
				numres = 2L,PACKAGE='pullword')
		OUT <- as.logical(OUT$numres)
	} else {
		if (combine) {
			OUT <- isUTF8(paste(string, collapse = ""))
		} else {
			OUT <- as.vector(sapply(string, isUTF8))
		}
	}
	return(OUT)
}

