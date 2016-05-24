#' Find Number of Unique Values
#'
#' Get the number of unique values in an input vector.
#'
#' @param vec A vector of any type.
#' @param na.rm If \code{TRUE}, remove missing values.
#'
#' @return Number of unique values.
#'
#' @export
#'
#' @examples
#' x <- c(1, 3, 1, 1, NA, 2, 2, 3, NA, NA, 1, 3, 1)
#' uniqlen(x)
#' uniqlen(x, na.rm=FALSE)

uniqlen <- function(vec, na.rm=TRUE) {
    if(na.rm && !is.null(vec)) {
       vec <- vec[!is.na(vec)]
    }
    length(unique(vec))
}
