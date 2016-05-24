#' Internal Ecospace Functions.
#'
#' Internal functions not intended to be called directly by users.
#'
#' Modified \code{unique} function that returns proper number of items for
#' filling of matrix when called within \code{apply}.
#'
#' @param x a vector or a data frame or an array or NULL.
#' @param length number of times (default=1) to repeat the dimensionally shorter
#'   state.
#' @param ... arguments for particular methods.
#'
#' @note When dealing with ecospace frameworks with multistate binary (numeric)
#'   character types and characters weighted by supplied species pools,
#'   sometimes all species will share the same state value for one of several
#'   states. (For example, perhaps all species are capable of sexual
#'   reproduction, but there is variation in whether some are exclusively sexual
#'   and some are hermaphroditic.) When this occurs, calling \code{apply} when
#'   choosing possible ecospace states will 'break the dimensionality' of the
#'   character matrix and convert it into a list. This function maintains matrix
#'   dimensionality, by repeating the dimensionally shorter unique state
#'   sufficient times to maintain equal length as found in other states.
#'
#' @seealso \code{\link[base]{unique}}
#'
#' @export
unique2 <- function(x, length=1, ...) {
  out <- unique(x, ...)
  if(length(out)==1) out <- rep(out, length)
  return(out)
  }
