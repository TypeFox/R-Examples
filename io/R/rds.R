#' @method qread rds
#' @export
qread.rds <- function(file, type, ...) {
	base::readRDS(file, ...)
}

#' @method qwrite rds
#' @export
qwrite.rds <- function(x, file, type, ...) {
	base::saveRDS(x, file, ...)
}
