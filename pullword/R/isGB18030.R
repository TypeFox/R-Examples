
##' Indicate whether the encoding of input string is GB18030.
##' 
##' @title Indicate whether the encoding of input string is GB18030.
##' @param string A character vector.
##' @param combine Whether to combine all the strings.
##' @return Logical value.
##' @author Jian Li <\email{rweibo@@sina.com}>

isGB18030 <- function(string, combine = FALSE)
{
	string <- .verifyChar(string)
	if (length(string)  == 1) {
		OUT <- .C("CWrapper_encoding_isgb18030", 
				characters = as.character(string),  
				numres = 2L,PACKAGE='pullword')
		OUT <- as.logical(OUT$numres)
	} else {
		if (combine) {
			OUT <- isGB18030(paste(string, collapse = ""))
		} else {
			OUT <- as.vector(sapply(string, isGB18030))
		}
	}
	return(OUT)
}

