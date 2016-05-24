#' @method qread data.frame
#' @export
qread.data.frame <- function(
	file, type, header=TRUE, sep=" ", check.names=TRUE, row.names=1, ...
) {
	read.table(file,
		header=header, sep=sep, check.names=check.names, row.names=row.names, ...)
}

#' @method qwrite data.frame
#' @export
qwrite.data.frame <- function(
	x, file, type, quote=TRUE, sep=" ", row.names=TRUE, col.names=NA, ...
) {
	write.table(x, file,
		quote=quote, sep=sep, row.names=row.names, col.names=col.names, ...)
}

#' @method qread dfm
#' @export
qread.dfm <- function(file, type, ...) {
	qread.data.frame(file)
}

#' @method qwrite dfm
#' @export
qwrite.dfm <- function(x, file, type, ...) {
	qwrite.data.frame(x, file)
}

