
check_string <- function(x) {
  stopifnot(is.character(x), length(x) == 1)
}

`%+%` <- function(lhs, rhs) {
  check_string(lhs)
  check_string(rhs)
  paste0(lhs, rhs)
}

chr <- as.character

int <- as.integer
