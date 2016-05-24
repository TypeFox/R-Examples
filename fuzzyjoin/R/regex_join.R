#' Join two tables based on a regular expression in one column
#' matching the other
#'
#' Join a table with a string column by a regular expression column
#' in another table
#'
#' @param x A tbl
#' @param y A tbl
#' @param by Columns by which to join the two tables
#' @param mode One of "inner", "left", "right", "full" "semi", or "anti"
#'
#' @seealso \code{\link{str_detect}}
#'
#' @examples
#'
#' library(dplyr)
#' library(ggplot2)
#' data(diamonds)
#'
#' diamonds <- tbl_df(diamonds)
#'
#' d <- data_frame(regex_name = c("^Idea", "mium", "Good"),
#'                 type = 1:3)
#'
#' # When they are inner_joined, only Good<->Good matches
#' diamonds %>%
#'   inner_join(d, by = c(cut = "regex_name"))
#'
#' # but we can regex match them
#' diamonds %>%
#'  regex_inner_join(d, by = c(cut = "regex_name"))
#'
#' @export
regex_join <- function(x, y, by = NULL, mode = "inner") {
  match_fun <- function(v1, v2) {
    stringr::str_detect(v1, v2)
  }

  fuzzy_join(x, y, by = by, match_fun = match_fun, mode = mode)
}


#' @rdname regex_join
#' @export
regex_inner_join <- function(x, y, by = NULL) {
  regex_join(x, y, by, mode = "inner")
}


#' @rdname regex_join
#' @export
regex_left_join <- function(x, y, by = NULL) {
  regex_join(x, y, by, mode = "left")
}


#' @rdname regex_join
#' @export
regex_right_join <- function(x, y, by = NULL) {
  regex_join(x, y, by, mode = "right")
}


#' @rdname regex_join
#' @export
regex_full_join <- function(x, y, by = NULL) {
  regex_join(x, y, by, mode = "full")
}


#' @rdname regex_join
#' @export
regex_semi_join <- function(x, y, by = NULL) {
  regex_join(x, y, by, mode = "semi")
}


#' @rdname regex_join
#' @export
regex_anti_join <- function(x, y, by = NULL) {
  regex_join(x, y, by, mode = "anti")
}
