fastdf <- function (list) {
  rows <- unique(unlist(lapply(list, NROW)))
  class(list) <- "data.frame"
  attr(list, "row.names") <- c(NA_integer_, -rows)
  list
}