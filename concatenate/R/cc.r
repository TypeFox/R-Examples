#' Comma Concatenation
#'
#' \code{cc} collapses text into a comma-separated list (the colloquial kind of
#' list). \code{cc_or} and \code{cc_and} insert "\code{or}" and "\code{or}"
#' before the last element.
#'
#' The \code{data.frame} method is dispatched when the first argument in
#' \code{...} is a \code{data.frame}. It operates row-wise. If there are
#' subsequent arguments to \code{cc} they are be ignored.
#' @param ... Character vectors or a \code{data.frame}.
#' @param oxford Whether to use the Oxford comma.
#' @return A length-one character vector in which each element in \code{...} is
#' separated by a comma (and a space).
#' @seealso \code{\link{cn}} for \code{cc} with (grammatical) number awareness
#' (like \code{ngettext}) and substitution (like \code{sprintf})
#' @examples
#' cc("hello", "world")
#'
#' a <- "one thing"
#' b <- "another"
#' cc_or(a, b)
#'
#' a <- "this"
#' b <- c("that", "the other")
#' cc_and(a, b)
#' @name cc
NULL

#' @rdname cc
#' @export
cc <- function(...) {
  s <- unlist(list(...))
  s <- trimws(s)
  paste(s, sep = ", ", collapse = ", ")
}

setGeneric("cc", signature = "...",
           function(...) standardGeneric("cc"))

#' @rdname cc
#' @export
setMethod("cc", "data.frame",
          function(...) {
            DF <- apply(list(...)[[1]], 1, cc)
            cc(DF)
          })

#' @rdname cc
#' @export
cc_or <- function(..., oxford = FALSE) {
  x = unlist(list(...))
  res <- cc(x[-length(x)])
  comma <- ifelse(isTRUE(oxford) && length(x) > 2, ",", "")
  or <- ifelse(length(x) > 1L, " or ", "")
  paste0(res, comma, or, x[length(x)])
}

#' @rdname cc
#' @export
cc_and <- function(..., oxford = FALSE) {
  x = unlist(list(...))
  res <- cc(x[-length(x)])
  comma <- ifelse(isTRUE(oxford) && length(x) > 2, ",", "")
  and <- ifelse(length(x) > 1L, " and ", "")
  paste0(res, comma, and, x[length(x)])
}

