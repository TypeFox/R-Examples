#' Find the least or most frequent occuring values in an atomic vector..
#'
#' frequent_values takes a vector and returns the n least or most occuring values.
#'
#' @param .x something.
#' @param n something.
#' @param type something.
#' @return  a vector.
#' @export
frequent_values <- function(.x, n = 5L, type = "most") {
  if(!is.atomic(.x)) stop(".x must be an atomic vector.")
  if(length(n) != 1L || !is.numeric(n) || any(n < 1)) {
    stop("n must be a single positive integer.")
  }
  if(!(type %in% c("most", "least"))) {
    stop("type takes only two values: 'most' or 'least'")
  }
  values <- sort(tapply(as.character(.x), factor(.x), length),
                 decreasing = TRUE)

  if (type == "most"){
    values_processed <-
      concat_rank_value_occurrence(r = unname(utils::head(rank(-values,
                                                               ties.method = "min",
                                                               na.last = FALSE), n)),
                                   n = utils::head(names(values), n),
                                   o = utils::head(unname(values), n))
  } else {
    values_processed <-
      concat_rank_value_occurrence(r = unname(rev(utils::tail(rank(-values,
                                                                   ties.method = "min",
                                                                   na.last = FALSE), n))),
                                   n = rev(utils::tail(names(values), n)),
                                   o = rev(utils::tail(unname(values), n)))
  }
  return(paste(values_processed, collapse = " "))
}
