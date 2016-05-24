#' ID By Row Number or Sequence Along
#' 
#' Generate a sequence of integers the 
#' \code{\link[base]{length}}/\code{\link[base]{ncol}} of an object.
#' 
#' @param x A dataframe, matrix, vector, or list object.
#' @param prefix A character string to use as a prefix. \code{FALSE} or 
#' \code{NULL} results in no prefix being used.  \code{TRUE} will utilize the 
#' prefix \code{"X."}.
#' @param pad logical.  If \code{TRUE} the beginning number will be padded with 
#' zeros.
#' @param \ldots Other arguments passed to \code{\link[qdapTools]{pad}}.
#' @return Returns a vector of sequential integers.
#' @keywords id
#' @export
#' @examples
#' id(list(1, 4, 6))
#' id(matrix(1:10, ncol=1))
#' id(mtcars)
#' id(mtcars, TRUE)
#' id("w")
#' id(mtcars, prefix="id-")
#' \dontrun{
#' library(qdap)
#' question_type(DATA.SPLIT$state, id(DATA.SPLIT, TRUE))
#' }
id <- function(x, prefix = FALSE, pad = TRUE, ...) {
  
    test1 <- dim(x)[1] > 1
    if (is.data.frame(x) | (!identical(logical(0), test1) && test1)) {
        ids <- seq_len(nrow(x))
    } else {
        ids <- seq_along(x)
    }
    if (pad) {
        ids <- pad(ids, ...)
    }
    if (!is.null(prefix) && is.character(prefix)) {
        ids <- paste0(prefix, ids)
    } 
    if (!is.null(prefix) && is.logical(prefix) && isTRUE(prefix)) {
        ids <- ds <- paste("X", ids, sep=".")
    } 	
    ids
}
