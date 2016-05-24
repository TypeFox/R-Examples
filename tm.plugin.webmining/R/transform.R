#' @title Enclose Text Content in HTML tags
#' @description Simple helper function which encloses text content of character
#' (or \code{\link[tm]{TextDocument}}) in HTML-tags. That way, HTML
#' content can be easier parsed by \code{\link[XML]{htmlTreeParse}}
#' @param x object of PlainTextDocument class
#' @export
#' @aliases encloseHTML.PlainTextDocument encloseHTML.character
encloseHTML <- function(x) UseMethod("encloseHTML", x)

#' @importFrom NLP content<-
#' @noRd
#' @export 
# FIXME: Could be done easier?? 
encloseHTML.PlainTextDocument <- function(x){
	content(x) <- sprintf("<html>%s</html>", x)
	x
} 

#' @title Remove non-ASCII characters from Text. 
#' @description This is a helper function to generate package data 
#' without non-ASCII character and omit the warning at R CMD check.
#' @param x object of PlainTextDocument class
#' @param fields specifies fields to be converted, defaults to fields = c("Content", "Heading", "Description")
#' @param from specifies encoding from which conversion should be done, defaults to "UTF-8"
#' @param to speciefies target encoding, defaults to "ASCII//TRANSLIT"
#' @export
#' @aliases removeNonASCII.PlainTextDocument
removeNonASCII <- function(x, fields = c("Content", "Heading", "Description"), from = "UTF-8", to = "ASCII//TRANSLIT")
	UseMethod("removeNonASCII", x)

#' @noRd
#' @export
removeNonASCII.PlainTextDocument <- function(x, fields = c("Content", "Heading", "Description"), from = "UTF-8", to = "ASCII//TRANSLIT"){
	if("Content" %in% fields){
		content(x) <- iconv(x, from, to)
	}
	for(fn in setdiff(fields, "Content")){
		meta(x, fn) <- iconv(meta(x, fn), from, to)
	}
	x
} 