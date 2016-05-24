#' @method qread xml
#' @export
qread.xml <- function(file, type, ...) {
	.check_package("XML");

	XML::xmlParse(file, ...)
}

#' @method qwrite xml
#' @export
qwrite.xml <- function(x, file, type, ...) {
	.check_package("XML");

	XML::saveXML(x, file, ...)
}
