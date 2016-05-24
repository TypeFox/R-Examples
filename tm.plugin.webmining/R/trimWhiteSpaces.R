#' @title Trim White Spaces from Text Document.
#' @description Transformation function, actually equal to stripWhiteSpace 
#' applicable for simple strings using Perl parser
#' @author Mario Annau
#' @param txt character
#' @seealso \code{\link{stripWhitespace}}
#' @export
trimWhiteSpaces <-
function(txt){
	txt <- sub("\\s+", "", txt, perl = TRUE)
	txt <- sub("\\s+$", "", txt, perl = TRUE)
	txt <- gsub("\\s\\s+", " ", txt, perl = TRUE)
	return(txt)
}

