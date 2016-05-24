#' fix_names
#'
#' Replace special characters in column names of a data.frame
#'
#' @param x a vector of column names
#' @param char substitute special char with this.
#'
#' @export
#'
#' @seealso make.names
#'
fix_names <- function (x, char = "."){

	y <- gsub("_", char, as.character(x))
	y <- gsub(" ", char, as.character(y))
	y <- gsub("\\.", char, as.character(y))

	return(y)
}
