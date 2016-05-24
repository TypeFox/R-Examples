#' @method qread list
#' @export
qread.list <- function(file, type, ...) {
	qread.yaml(file, type, ...)
}

#' @method qwrite list
#' @export
qwrite.list <- function(x, file, type, ...) {
	qwrite.yaml(x, file, type, ...)
}

#' @method qread lst
#' @export
qread.lst <- qread.list;

#' @method qwrite lst
#' @export
qwrite.lst <- qwrite.list;

