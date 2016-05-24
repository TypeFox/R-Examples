#' @rdname bpa
#' @export
get_pattern <- function(x, show_ws = TRUE, ws_char = "w") {
  if (!is.character(x)) {
    x <- as.character(x)
  }
  x <- gsub("[a-z]", "a", x)
  x <- gsub("[A-Z]", "A", x)
  x <- gsub("[0-9]", "9", x)
  if (show_ws) {
    x <- gsub("\\s", ws_char, x)
  }
  x
}

#' Basic Pattern Analysis
#'
#' Perform a basic pattern analysis
#'
#' @param x A data frame or character vector.
#' @param unique_only Logical indicating whether or not to only show the unique
#'   patterns. Default is \code{TRUE}.
#' @param show_ws Logical indicating whether or not to show whitespace
#'   using a special character. Default is \code{TRUE}.
#' @param ws_char Character string to use to depict whitespace when 
#'   \code{show_ws = TRUE}.
#' @param useNA Logical indicating whether to include \code{NA} values in the
#'   table. See \code{\link{table}} for details.
#' @param ... Additional optional arguments to be passed onto \code{llply}.
#' @rdname bpa
#' @export
#' @examples
#' basic_pattern_analysis(iris)
#' basic_pattern_analysis(iris, unique_only = TRUE)
basic_pattern_analysis <- function(x, unique_only = FALSE,
                                   show_ws = TRUE, ws_char = "w",
                                   useNA = c("no", "ifany", "always"), ...) {
  UseMethod("basic_pattern_analysis")
}


#' @rdname bpa
#' @export
basic_pattern_analysis.default <- function(x, unique_only = FALSE,
                                           show_ws = TRUE, ws_char = "w",
                                           useNA = c("no", "ifany", "always"),
                                           ...) {
  useNA <- match.arg(useNA)
  if (unique_only) {
    table(get_pattern(x, show_ws = show_ws, ws_char = ws_char), useNA = useNA)
  } else {
    get_pattern(x, show_ws = show_ws, ws_char = ws_char)
  }
}


#' @rdname bpa
#' @importFrom plyr llply
#' @export
basic_pattern_analysis.data.frame <- function(x, unique_only = FALSE,
                                              show_ws = TRUE, ws_char = "w",
                                              useNA = c("no", "ifany", "always"),
                                              ...) {
  useNA <- match.arg(useNA)
  z <- llply(x, basic_pattern_analysis.default, unique_only = unique_only,
             show_ws = show_ws, ws_char = ws_char, useNA = useNA, ...)
  if (unique_only) {
    z
  } else {
    data.frame(z)
  }
}


#' @rdname bpa
#' @export
bpa <- function(x, ...){
  basic_pattern_analysis(x, ...)
}
