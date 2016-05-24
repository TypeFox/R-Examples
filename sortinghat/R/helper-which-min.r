#' Helper function that determines which element in a vector is the minimum. Ties
#' can be broken randomly or via first/last ordering.
#'
#' The \code{which_min} function is intended to be an alternative to the base
#' \code{which.min} function when a specific tie-breaking method is necessary.
#'
#' @export
#' @param x vector
#' @param break_ties method to break ties. The \code{random} method selects the
#' index of the minimum elements randomly, while the \code{first} and \code{last}
#' options imply that the first or last instance of the minimum element will be
#' chosen, respectively.
#' @return location of the minimum element in the vector \code{x}. If there is a
#' tie, we break the tie with the method specified in \code{break_ties}.
#' @examples
#' set.seed(42)
#' z <- runif(5)
#' z <- c(z[1], z[1], z)
#'
#' which_min(z)
#' which_min(z, break_ties = "first")
#' which_min(z, break_ties = "last")
which_min <- function(x, break_ties = c("random", "first", "last")) {
  break_ties <- match.arg(break_ties)
  min_x <- min(x)
  which_min_x <- sort(which(x == min_x))

  # If there is a tie, we break the tie according to the specified method.
  if (length(which_min_x) > 1) {
    if (break_ties == "random") {
      which_min_x <- sample(which_min_x, size = 1)
    } else if (break_ties == "first") {
      which_min_x <- head(which_min_x, 1)
    } else if (break_ties == "last") {
      which_min_x <- tail(which_min_x, 1)
    }
  }
  which_min_x
}

