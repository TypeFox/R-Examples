
#' Concentration indices
#'
#' This is a proof-of-concept implementation of concentration measures
#' for data with an exogenous order.
#'
#' @param x Data vector
#'
#' @return The concenctration index.
#'
#' @examples
#'   ## See demo("concentration-bundesliga") for a real-world
#'   ## application.
#'
#'   set.seed(1234)
#'   x <- sample(100000, 10)
#'
#'   ## Classical concentration indices:
#'   herfindahl(x)
#'   rosenbluth(x)
#'
#'   ## Exogenous order is available:
#'   o <- sample(10)
#'   exogeny(x, o)
#'
#' @family concentration
#' @rdname concentration_indices
#' @aliases concentration_indices
#'
#' @export
herfindahl <- function(x) {
  y <- sumex(x, ex = order(x))^2
  y <- sum(y)

  y
}


#' @rdname concentration_indices
#' @export
rosenbluth <- function(x) {
  n <- length(x)

  y <- sumex(x, ex = order(x, decreasing = TRUE))
  y <- sum((1:n) * y)
  y <- 1 / ((2 * y) - 1)

  y
}


#' @param ex Order of the data vector
#' @rdname concentration_indices
#' @export
exogeny <- function(x, ex = order(x)) {
  n <- length(x)

  y <- sumex(x, ex = ex)
  y <- sum((1:n) * y)
  y <- 1 / ((2 * y) - 1)

  y
}


#' Concentration ratios
#'
#' This is a proof-of-concept implementation of concentration ratios
#' for data with an exogenous order.
#'
#' @param x Data vector
#' @param g Number of "first" data (according to the order)
#' @param ex Exogenous order of the data vector
#' @return The concentration ratio value
#'
#' @examples
#'   set.seed(1234)
#'   x <- sample(100000, 10)  # Data
#'   o <- sample(10)          # Exogenous order
#'
#'   concentration_ratio(x, 3)
#'   concentration_ratio(x, 3, ex = o)
#'
#' @family concentration
#'
#' @export
concentration_ratio <- function(x, g, ex = order(x, decreasing = TRUE)) {
  y <- sumex(x, ex = ex)
  y <- sum(y[1:g])

  y
}


sumex <- function(x, ex) {
  x <- x[ex]
  x / sum(x)
}


