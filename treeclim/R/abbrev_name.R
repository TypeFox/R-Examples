##' Shorten month names
##' 
##' Shorten month names that are longer than 43 characters by joining
##' the first and last 20 characters with "...".
##' @param x a vector of month names as returned by format_month()
##' @return a vector of month names shortened to 43 characters
##' @keywords manip, internal
abbrev_name <- function(x) {
  .x <- x
  n <- length(x)
  for (i in 1:n) {
    k <- nchar(x[i])
    if (k > 43) {
      ## take first and last 20 chars
      pref <- substr(x[i], 1, 20)
      suff <- substr(x[i], k-19, k)
      .x[i] <- paste(pref, suff, sep = "...")
      .x[i] <- gsub("\\.\\.\\.\\.", "...", .x[i])
    }
  }
  .x
}
