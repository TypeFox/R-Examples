#' Determine input-output support for data or file type
#'
#' This function returns whether a type is supported by
#' \code{\link{qread}} or \code{\link{qwrite}}.
#'
#' @param type  data or file type
#' @return a \code{data.frame} with logical entries;
#'         \code{TRUE} if type is supported, \code{FALSE} otherwise
#' @export
#'
#' @examples
#' io_supported("rds")
#'
io_supported <- function(type) {
	type <- tolower(type);
	data.frame(
		qread = type %in% sub("^qread\\.", "", as.character(methods(qread))),
		qwrite = type %in% sub("^qwrite\\.", "", as.character(methods(qwrite))),
		row.names=type
	)
}
