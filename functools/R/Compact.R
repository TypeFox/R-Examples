#' Filter NA and NULL values out of a vector, list, or data.frame.
#'
#' \code{Compact()} takes a vector \code{.x} and returns it with all NULL and
#' NA values filtered out.
#'
#' @param .x A vector.
#' @return Vector .x but with all NULL and NA values filtered out.
#' @examples
#' # Removes all null elements from a vector:
#' a <- list(NULL, 1, 5, NULL)
#' Compact(a)
#'
#' b <- c(1, 2, 0, 4, NULL, 1, 3, NULL)
#' Compact(b)
#'
#' @export
Compact <- function(.x) return(Filter(Existy, .x))
