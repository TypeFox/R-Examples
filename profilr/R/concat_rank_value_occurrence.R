#' Wrapper for concatenating rank, value, and number of occurrences.
#'
#'concat_rank_value_occurrence takes vectors representing rank, value, and number of occurrences and returns them in a concatenated string. This is just a wrapper around paste0().
#'
#' @param r something.
#' @param n something.
#' @param o something.
#' @return  a vector.
#' @examples
#' # Example
#' concat_rank_value_occurrence(1, "a", 100)
#' @export
concat_rank_value_occurrence <- function(r, n, o) {
  paste0("(", r, ")", " ", n, " [", o, "]", "\n", sep = "")
}
