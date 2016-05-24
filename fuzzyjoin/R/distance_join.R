#' Join two tables based on a distance metric of one or more columns
#'
#' This differs from \code{\link{difference_join}} in that it considers
#' all of the columns together when computing distance. This allows it
#' to use metrics such as Euclidean or Manhattan that depend on multiple
#' columns. Note that if you are computing with longitude or latitude,
#' you probably want to use \code{\link{geo_join}}.
#'
#' @param x A tbl
#' @param y A tbl
#' @param by Columns by which to join the two tables
#' @param max_dist Maximum distance to use for joining
#' @param method Method to use for computing distance, either euclidean (default)
#' or manhattan.
#' @param mode One of "inner", "left", "right", "full" "semi", or "anti"
#'
#' @examples
#'
#' library(dplyr)
#'
#' head(iris)
#' sepal_lengths <- data_frame(Sepal.Length = c(5, 6, 7),
#'                             Sepal.Width = 1:3)
#'
#' iris %>%
#'   distance_inner_join(sepal_lengths, max_dist = 2)
#'
#' @export
distance_join <- function(x, y, by = NULL, max_dist = 1,
                          method = c("euclidean", "manhattan"),
                          mode = "inner") {
  method <- match.arg(method)

  match_fun <- function(v1, v2) {
    if (method == "euclidean") {
      d <- sqrt(rowSums((v1 - v2) ^ 2))
    } else if (method == "manhattan") {
      d <- rowSums(abs(v1 - v2))
    }
    d <= max_dist
  }

  fuzzy_join(x, y, multi_by = by, multi_match_fun = match_fun, mode = mode)
}


#' @rdname distance_join
#' @export
distance_inner_join <- function(x, y, by = NULL, method = "euclidean", max_dist = 1) {
  distance_join(x, y, by, max_dist = max_dist, method = method, mode = "inner")
}


#' @rdname distance_join
#' @export
distance_left_join <- function(x, y, by = NULL, method = "euclidean", max_dist = 1) {
  distance_join(x, y, by, max_dist = max_dist, method = method, mode = "left")
}


#' @rdname distance_join
#' @export
distance_right_join <- function(x, y, by = NULL, method = "euclidean", max_dist = 1) {
  distance_join(x, y, by, max_dist = max_dist, method = method, mode = "right")
}


#' @rdname distance_join
#' @export
distance_full_join <- function(x, y, by = NULL, method = "euclidean", max_dist = 1) {
  distance_join(x, y, by, max_dist = max_dist, method = method, mode = "full")
}


#' @rdname distance_join
#' @export
distance_semi_join <- function(x, y, by = NULL, method = "euclidean", max_dist = 1) {
  distance_join(x, y, by, max_dist = max_dist, method = method, mode = "semi")
}


#' @rdname distance_join
#' @export
distance_anti_join <- function(x, y, by = NULL, method = "euclidean", max_dist = 1) {
  distance_join(x, y, by, max_dist = max_dist, method = method, mode = "anti")
}
