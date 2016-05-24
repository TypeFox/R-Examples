#' @title Wrapper/Convenience function to ensure right encoding for different Platforms
#' @description Depending on specified type one of the following parser functions is called:
#' \describe{
#' \item{XML}{\code{\link{xmlInternalTreeParse}}}
#' \item{HTML}{\code{\link{htmlTreeParse}}}
#' \item{JSON}{\code{\link{fromJSON}}}
#' }
#' @param ... arguments to be passed to specified parser function
#' @param asText defines if input should be treated as text/character, default to TRUE
#' @param type either "XML", "HTML" or "JSON". Defaults to "XML"
#' @export
parse <- function(..., asText = TRUE, type = c("XML", "HTML", "JSON")){
	parsetype <- match.arg(type)
	encoding <- switch(.Platform$OS.type,
						unix = "UTF-8",
						windows = "latin1")
	parser <- switch(parsetype,
						XML = xmlInternalTreeParse,
						HTML = htmlTreeParse,
						JSON = fromJSON)
	parser(..., encoding = encoding, asText = asText)
}
