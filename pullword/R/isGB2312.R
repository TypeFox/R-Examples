
##' Indicate whether the encoding of input string is GB2312.
##' 
##' @title Indicate whether the encoding of input string is GB2312.
##' @param string A character vector.
##' @param combine Whether to combine all the strings.
##' @return Logical value.
##' @author Jian Li <\email{rweibo@@sina.com}>

isGB2312 <- function(string, combine = FALSE)
{
	string <- .verifyChar(string)
	if (length(string)  == 1) {
		OUT <- .C("CWrapper_encoding_isgb2312", 
				characters = as.character(string),  
				numres = 2L,PACKAGE='pullword')
		OUT <- as.logical(OUT$numres)
	} else {
		if (combine) {
			OUT <- isGB2312(paste(string, collapse = ""))
		} else {
			OUT <- as.vector(sapply(string, isGB2312))
		}
	}
	return(OUT)
}

