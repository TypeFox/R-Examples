#' setNames for arrays...
#'
#' @param array.object Array
#' @param list.of.dimension.names List of character vectors with which to name each dimension of the array
#' @return Named array
#' @examples
#' setDimNames(matrix(1:4,2,2), list(c("Cat", "Dog"), c("Name", "Weight")))
#' @export
setDimNames <- function(array.object, list.of.dimension.names) {
	dimnames(array.object) <- list.of.dimension.names
	array.object
}

#' Capitalise words in character vector
#'
#' @param x Character vector 
#' @return Character vector
#' @examples
#' simpleCap(c("a simple test", "Another-test"))
#' @export
simpleCap <- function(x) {
	s <- strsplit(x, " ")
	sapply(
		s,
		function(x) paste(
			toupper(substring(x, 1,1)), 
			substring(x, 2),
			sep="", collapse=" "
		)
	)
}

