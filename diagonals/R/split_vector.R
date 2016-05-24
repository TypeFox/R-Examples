#'Split Vector
#'
#' @name split_vector
#' @param x a numeric or character vector
#' @param steps the number of steps
#' @param size the size of the step
#' @param replacement value to be inserted on the diagonal, by default this is
#'   zero (0).
#' @details Either steps or size is expected to be provided.
#' @export


split_vector <- function( x, steps = NULL, size = NULL, replacement = 0 ) {

  # determine the size of the step
  if ( !is.null(size) ) {

    # coerce to integer
    size <- as.integer(size)

  } else if ( is.null(size) & !is.null(steps) ) {

    # calculate size
    size <- as.integer( length(x) %/% steps )

  } else if (is.null(steps) & is.null(size) ) {

    # issue warning
    warning(paste("Both steps and size parameters are NULL, setting step size to 1 (one). ") )

    # set to unit
    size <- 1L

  }

  split( x, ceiling( seq_along(x)/size ) )

}
