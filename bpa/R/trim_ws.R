#' Remove Leading/Trailing Whitespace
#'
#' Remove leading and/or trailing whitespace from character strings.
#'
#' @param x A data frame or vector.
#' @param which A character string specifying whether to remove both leading and
#'   trailing whitespace (default), or only leading (\code{"left"}) or trailing
#'   (\code{"right"}). Can be abbreviated.
#' @export
#' @examples
#' # Toy example
#' d <- data.frame(x = c(" a ", "b ", "c"),
#'                 y = c("   1 ", "2", " 3"),
#'                 z = c(4, 5, 6))
#' print(d)  # print data as is
#' trim_ws(d)  # print data with whitespace trimmed off
#' sapply(trim_ws(d), class)  # check that column types are preserved
trim_ws <- function(x, which = c("both", "left", "right")) {
  UseMethod("trim_ws")
}


#' @export
trim_ws.default <- function(x, which = c("both", "left", "right")) {
  # trimws(x, ...)
  .which <- match.arg(which)
  if (.which == "both") {
    gsub("(^\\s+)|(\\s+$)", "", x)
  } else if (.which == "left") {
    gsub("(^\\s+)", "", x)
  } else {
    gsub("(\\s+$)", "", x)
  }
}


#' @export
trim_ws.data.frame <- function(x, which = c("both", "left", "right")) {
  .which <- match.arg(which)
  as.data.frame(lapply(x, function(y) {
    if (is.numeric(y)) y else trim_ws(y, which = .which)
  }), stringsAsFactors = FALSE)
}
