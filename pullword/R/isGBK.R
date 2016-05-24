
##' Indicate whether the encoding of input string is GBK.
##' 
##' @title Indicate whether the encoding of input string is GBK.
##' @param string A character vector.
##' @param combine Whether to combine all the strings.
##' @return Logical value.
##' @author Jian Li <\email{rweibo@@sina.com}>

isGBK <- function(string, combine = FALSE)
{
	string <- .verifyChar(string)
	if (length(string)  == 1) {
		OUT <- .C("CWrapper_encoding_isgbk", 
				characters = as.character(string),  
				numres = 2L,PACKAGE='pullword')
		OUT <- as.logical(OUT$numres)
	} else {
		if (combine) {
			OUT <- isGBK(paste(string, collapse = ""))
		} else {
			OUT <- as.vector(sapply(string, isGBK))
		}
	}
	return(OUT)
}

