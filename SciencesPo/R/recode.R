#' @encoding UTF-8
#' @title Recode or Replace Values With New Values
#'
#' @description Recodes a value or a vector of values.
#'
#' @param x The vector whose values will be recoded.
#' @param from a vector of the items to recode.
#' @param to a vector of replacement values.
#' @param warn A logical to print a message if any of the
#'  old values are not actually present in \code{x}.
#' @keywords Manipulation
#'
#' @examples
#' x <- LETTERS[1:5]
#' recode(x, c("B", "D"), c("Beta", "Delta"))
#'
#' # On numeric vectors
#' x <- c(1, 4, 5, 9)
#' recode(x, from = c(1, 4, 5, 9), to = c(10, 40, 50, 90))
#'
#' @export
`recode` <- function(x, from, to, warn=TRUE) UseMethod("recode")
NULL


#' @rdname recode
#' @export
`recode.default` <- function(x, from, to, warn = TRUE) {
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x)) {
    stop("`x` must be an atomic vector.")
  }

  if (is.factor(x)) {
    # If x is a factor, call self but operate on the levels
    levels(x) <- recode(levels(x), from, to, warn)
    return(x)
  }

  map_x <- match(x, from)
  map_x_NA  <- is.na(map_x)

  # index of items in `from` that were found in `x`
  from_unique <- sort(unique(map_x))
  if (warn && length(from_unique) != length(from)) {
    message("The following `from` values were not present in `x`: ",
     paste(from[!(1:length(from) %in% from_unique) ], collapse = ", "))
  }
  x[!map_x_NA] <- to[map_x[!map_x_NA]]
  return(x)
}
NULL


