#' Get Names of Attributes
#'
#' Get the names of the attributes of an input object.
#'
#' @param obj Any object.
#'
#' @return Vector of character strings with the names of the attributes.
#'
#' @export
#'
#' @examples
#' x <- matrix(1:100, ncol=5)
#' colnames(x) <- LETTERS[1:5]
#' attrnames(x)

attrnames <- function(obj) {
  names(attributes(obj))
}
