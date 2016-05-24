#' Reformat time.
#'
#' Functions perform interconversion between "HH:MM:SS" format and seconds.
#'
#' @param x either a character string of the form "HH:MM:SS" ("HH" is optional)
#'   or numeric \strong{seconds} values.
#'
#' @return seconds value(s) for \emph{from}, and "HH:MM:SS" character string(s)
#'   for \emph{to}.
#'
#' @examples
#' x <- c("00:21:05", "25:51", NA, "00:26:01.1", "01:05:02.0")
#' x <- convert_from_time(x)
#' print(x)
#' x <- convert_to_time(x)
#' print(x)
#'
#' @name convert_time
NULL
# ------------------------------------------------------------------------------
#' @rdname convert_time
#' @export
convert_from_time <- function(x) {
  stopifnot(is.character(x))
  x <- suppressWarnings(strsplit(x, ":"))
  x <- vapply(x, FUN.VALUE = numeric(1), function(i) {
    sum(as.numeric(rev(i)) * c(1, 60, 60 ^ 2)[seq_len(length(i))])
  })
  x
}
#' @rdname convert_time
#' @export
convert_to_time <- function(x) {
  stopifnot(is.numeric(x))
  h  <- x / (60 ^ 2)
  m  <- (h - (h <- floor(h))) * 60
  s  <- (m - (m <- floor(m))) * 60
  o  <- sprintf("%02.0f:%02.0f:%02.1f", h, m, s)
  NA -> o[grep("NA", o, ignore.case = TRUE)]
  o
}
