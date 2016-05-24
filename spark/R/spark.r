
spark_ticks <- c("\u2581", "\u2582", "\u2583", "\u2584",
                 "\u2585", "\u2586", "\u2587", "\u2588")

spark_ticks_fallback <- c("_", "~", "^")

#' Generic spark line method for matrices or data frames, etc.
#'
#' @param data The data to plot.
#' @param width The width (number of characters) of the output.
#'   \sQuote{data} means that it will match the length of the data.
#'   \sQuote{screen} means that it will be scaled to match the
#'   width of the screen. \sQuote{auto} means \sQuote{data}
#'   if the length of the data is not longer than the screen width,
#'   and \sQuote{screen} otherwise.
#' @param print Whether to show the result on the screen. (If \code{FALSE},
#'   it will be only returned.)
#' @param ... Additional arguments to the various methods.
#' @return Characacter vector containing the spark line(s).
#' @family spark
#' @export

spark <- function(data, width = c("data", "auto", "screen"), print = TRUE,
                  ...)
  UseMethod("spark")


is_utf8 <- function() {

  Sys.getlocale("LC_CTYPE") %>%
    grepl(pattern = "UTF-8")

}


.onLoad <- function(libname, pkgname) {

  if (!is_utf8()) {
    assign("spark_ticks", spark_ticks_fallback, envir = asNamespace(pkgname))
  }

}
