#' Join two tables based on absolute difference between their columns
#'
#' @param x A tbl
#' @param y A tbl
#' @param by Columns by which to join the two tables
#' @param max_dist Maximum distance to use for joining
#' @param mode One of "inner", "left", "right", "full" "semi", or "anti"
#'
#' @examples
#'
#' library(dplyr)
#'
#' head(iris)
#' sepal_lengths <- data_frame(Sepal.Length = c(5, 6, 7), Type = 1:3)
#'
#' iris %>%
#'   difference_inner_join(sepal_lengths, max_dist = .5)
#'
#' @export
difference_join <- function(x, y, by = NULL, max_dist = 1, mode = "inner") {
  match_fun <- function(v1, v2) {
    abs(v1 - v2) <= max_dist
  }

  fuzzy_join(x, y, by = by, match_fun = match_fun, mode = mode)
}


#' @rdname difference_join
#' @export
difference_inner_join <- function(x, y, by = NULL, max_dist = 1) {
  difference_join(x, y, by, max_dist = max_dist, mode = "inner")
}


#' @rdname difference_join
#' @export
difference_left_join <- function(x, y, by = NULL, max_dist = 1) {
  difference_join(x, y, by, max_dist = max_dist, mode = "left")
}


#' @rdname difference_join
#' @export
difference_right_join <- function(x, y, by = NULL, max_dist = 1) {
  difference_join(x, y, by, max_dist = max_dist, mode = "right")
}


#' @rdname difference_join
#' @export
difference_full_join <- function(x, y, by = NULL, max_dist = 1) {
  difference_join(x, y, by, max_dist = max_dist, mode = "full")
}


#' @rdname difference_join
#' @export
difference_semi_join <- function(x, y, by = NULL, max_dist = 1) {
  difference_join(x, y, by, max_dist = max_dist, mode = "semi")
}


#' @rdname difference_join
#' @export
difference_anti_join <- function(x, y, by = NULL, max_dist = 1) {
  difference_join(x, y, by, max_dist = max_dist, mode = "anti")
}
