#' @method qread json
#' @export
qread.json <- function(file, type, ...) {
	.check_package("jsonlite");

	jsonlite::unserializeJSON(paste(readLines(file), collapse="\n"), ...)
}

#' @method qwrite json
#' @export
qwrite.json <- function(x, file, type, append=FALSE, ...) {
	.check_package("jsonlite");

	cat(jsonlite::serializeJSON(x, ...), sep="", file=file, append=append)
}
