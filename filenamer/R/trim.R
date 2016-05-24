#' Trim extensions from a file name
#'
#' This function trims extensions from a file name.
#' 
#' @param x     a \code{character} or a \code{filename}
#' @param n     number of extensions to trim off the end
#' @return modified object of the original type
#' @export
#' 
#' @examples
#' x <- "path/data.txt.gz"
#' print(trim_ext(x))
#'
trim_ext <- function(x, n) UseMethod("trim_ext");

#' @export
trim_ext.filename <- function(x, n=1) {
	if (length(x$ext) <= n) {
		x$ext <- NULL;
	} else {
		x$ext <- x$ext[1:(length(x$ext)-n)];
	}
	x
}

#' @export
trim_ext.character <- function(x, n=1) {
	as.character(trim_ext.filename(as.filename(x), n))
}
