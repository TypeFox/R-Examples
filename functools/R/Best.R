#' Find the best value in a vector.
#'
#' \code{Best()} takes a vector \code{.x} and a binary predicate function
#' \code{.f} and returns the result of \code{.f} reduced over \code{.x}.
#'
#' @param .x A vector.
#' @param .f A binary predicate function.
#' @return The best value in that vector, as determined by the binary predicate function.
#' @family aggregate functionals
#' @examples
#' # Simulate the behavior of max with numerics
#' Best(1:10, function(x, y) return(x > y))
#' # Simulate the behavior of min with numerics
#' Best(1:10, function(x, y) return(x < y))
#' # This comparison function prefers values that begin with l
#' Best(letters, function(x, y) return(x[1] == "l"))
#' @export
Best <- function(.x, .f) {
  .f <- match.fun(.f)
  return(Reduce(function(.y, .z) ifelse(.f(.y, .z), .y, .z), .x))
}
