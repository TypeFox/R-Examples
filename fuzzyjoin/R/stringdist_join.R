#' Join two tables based on fuzzy string matching of their columns
#'
#' Join two tables based on fuzzy string matching of their columns.
#' This is useful, for example, in matching free-form inputs in
#' a survey or online form, where it can catch misspellings and
#' small personal changes.
#'
#' @param x A tbl
#' @param y A tbl
#' @param by Columns by which to join the two tables
#' @param max_dist Maximum distance to use for joining
#' @param ignore_case Whether to be case insensitive (default yes)
#' @param method Method for computing string distance, see
#' \code{stringdist-methods} in the stringdist package.
#' @param mode One of "inner", "left", "right", "full" "semi", or "anti"
#' @param ... Arguments passed on to \code{\link{stringdist}}
#'
#' @details If \code{method = "soundex"}, the \code{max_dist} is
#' automatically set to 0.5, since soundex returns either a 0 (match)
#' or a 1 (no match).
#'
#' @examples
#'
#' library(dplyr)
#' library(ggplot2)
#' data(diamonds)
#'
#' d <- data_frame(approximate_name = c("Idea", "Premiums", "Premioom",
#'                                      "VeryGood", "VeryGood", "Faiir"),
#'                 type = 1:6)
#'
#' # no matches when they are inner-joined:
#' diamonds %>%
#'   inner_join(d, by = c(cut = "approximate_name"))
#'
#' # but we can match when they're fuzzy joined
#' diamonds %>%
#'  stringdist_inner_join(d, by = c(cut = "approximate_name"))
#'
#' @export
stringdist_join <- function(x, y, by = NULL, max_dist = 2,
                            method = c("osa", "lv", "dl", "hamming", "lcs", "qgram",
                                       "cosine", "jaccard", "jw", "soundex"),
                            mode = "inner",
                            ignore_case = FALSE, ...) {
  method <- match.arg(method)

  if (method == "soundex") {
    # soundex always returns 0 or 1, so any other max_dist would
    # lead either to always matching or never matching
    max_dist <- .5
  }

  match_fun <- function(v1, v2) {
    if (ignore_case) {
      v1 <- stringr::str_to_lower(v1)
      v2 <- stringr::str_to_lower(v2)
    }

    # shortcut for Levenshtein-like methods: if the difference in
    # string length is greater than the maximum string distance, the
    # edit distance must be at least that large

    # length is much faster to compute than string distance
    if (method %in% c("osa", "lv", "dl")) {
      length_diff <- abs(stringr::str_length(v1) - stringr::str_length(v2))
      include <- length_diff <= max_dist

      ret <- logical(length(v1))

      dists <- stringdist::stringdist(v1[include], v2[include], method = method, ...)
      ret[include] <- dists <= max_dist
      ret
    } else {
      # have to compute them all
      dists <- stringdist::stringdist(v1, v2, method = method, ...)
      dists <= max_dist
    }
  }

  fuzzy_join(x, y, by = by, mode = mode, match_fun = match_fun)
}


#' @rdname stringdist_join
#' @export
stringdist_inner_join <- function(x, y, by = NULL, ...) {
  stringdist_join(x, y, by, mode = "inner", ...)
}


#' @rdname stringdist_join
#' @export
stringdist_left_join <- function(x, y, by = NULL, ...) {
  stringdist_join(x, y, by, mode = "left", ...)
}


#' @rdname stringdist_join
#' @export
stringdist_right_join <- function(x, y, by = NULL, ...) {
  stringdist_join(x, y, by, mode = "right", ...)
}


#' @rdname stringdist_join
#' @export
stringdist_full_join <- function(x, y, by = NULL, ...) {
  stringdist_join(x, y, by, mode = "full", ...)
}


#' @rdname stringdist_join
#' @export
stringdist_semi_join <- function(x, y, by = NULL, ...) {
  stringdist_join(x, y, by, mode = "semi", ...)
}


#' @rdname stringdist_join
#' @export
stringdist_anti_join <- function(x, y, by = NULL, ...) {
  stringdist_join(x, y, by, mode = "anti", ...)
}
