#' Insert tag or extension and coerce to character
#'
#' This fucntion inserts a tag or extension into a file name and returns
#' a \code{charcter} vector.
#'
#' @param x    a \code{filename} or \code{character}
#' @param ...  arguments passed to \code{\link{insert}}
#' @return a \code{character} vector
#' @export
#'
#' @examples
#' x <- "data.txt"
#' y <- tag(x, "qc")
#' print(y)
#' f <- as.filename(x)
#' g <- tag(f, "qc")
#' print(g)
#'
tag <- function(x, ...) UseMethod("tag");

#' @export
tag.filename <- function(x, ...) {
	as.character( insert.filename(x, ...) )
}

#' @export
tag.character <- function(x, ...) {
	insert.character(x, ...)
}
