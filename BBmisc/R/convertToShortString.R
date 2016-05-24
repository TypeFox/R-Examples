#' @title Converts any R object to a descriptive string so it can be used in messages.
#'
#' @description
#' Atomics: If of length 0 or 1, they are basically printed as they are.
#' Numerics are formated with \code{num.format}.
#' If of length greater than 1, they are collapsed witd \dQuote{,} and clipped.
#' so they do not become excessively long.
#'
#' All others: Currently, only their class is simply printed
#' like \dQuote{<data.frame>}.
#'
#' Lists: The mechanism above is applied (non-recursively) to their elements.
#' The result looks like this:
#' \dQuote{a = 1, <unamed> = 2, b = <data.frame>, c = <list>}.
#'
#' @param x [any]\cr
#'   The object.
#' @param num.format [\code{character(1)}]\cr
#'   Used to format numerical scalars via \code{\link{sprintf}}.
#'   Default is \dQuote{\%.4g}.
#' @param clip.len [\code{integer(1)}]\cr
#'   Used clip atomic vectors via \code{\link{clipString}}.
#'   Default is 15.
#' @return [\code{character(1)}].
#' @export
#' @examples
#' convertToShortString(list(a = 1, b = NULL, "foo", c = 1:10))
convertToShortString = function(x, num.format = "%.4g", clip.len = 15L) {

  # convert non-list object to string
  convObj = function(x) {
    if (is.null(x))
      return("NULL")
    if (is.atomic(x)) {
      if (length(x) == 0L) {
        sprintf("%s(0)", getClass1(x))
      } else if (length(x) == 1L) {
        if (is.double(x))
          sprintf(num.format, x)
        else
          collapse(x)
      } else {
        clipString(collapse(sapply(x, convertToShortString), ","), clip.len)
      }
    } else {
      paste("<", getClass1(x), ">", sep = "")
    }
  }

  # handle only lists and not any derived data types
  if (getClass1(x) == "list") {
    if (length(x) == 0L)
      return("")
    ns = names2(x, missing.val = "<unnamed>")
    ss = lapply(x, convObj)
    collapse(paste(ns, "=", ss, sep = ""), ", ")
  } else {
    convObj(x)
  }
}
