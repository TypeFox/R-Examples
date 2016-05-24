if(getRversion() >= "2.15.1") globalVariables(c(".data", "copy", ".N"))
#' @encoding UTF-8
#' @title Number of the Current Observation
#'
#' @description Computes the number of the current observation. It's provides similar functionality of the Stata "_n" built-in system variable.
#' @param .data The input \code{data.frame}.
#' @param by The variables that should be treated as ID variables. Defaults to \code{NULL}, at which point none variable is  used as ID variable.
#' @examples
#' df <- data.frame(A = c("a", "a", "a", "b", "b"),
#'                  B = c(1, 1, 1, 1, 1), values = 1:5);
#' df;
#'
#' .n(df);
#'
#' .n(df, by = "A");
#'
#' # Numeric also works fine:
#' .n(df, by = 1);
#'
#' @importFrom data.table := is.data.table as.data.table
#' @export
`.n` <- function(.data, by = NULL) {
  if (!is.data.table(.data)) .data <- as.data.table(.data)
  else .data <- copy(.data)

  if (is.numeric(by)) by <- names(.data)[by]
  if (is.null(by)) by <- names(.data)

  `.n` <- .N  <- NULL

  if (any(duplicated(.data, by = by))) {
    .data[, `.n` := sequence(.N), by = by]
  } else {
    .data[, `.n` := sequence(.N)]
  }
}
NULL
