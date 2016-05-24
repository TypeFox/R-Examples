#' Number-aware Strings with Substitution
#'
#' \code{cn} combines grammatical number awareness as in \code{\link{ngettext}}
#' with \code{\link{sprintf}}-like substitution for comma-concatenated text.
#'
#' Like \code{ngettext}, this function returns one string to be used with a
#' singular referent and another with a plural referent. \code{cn} chooses
#' between the two based on the length of its first argument, \code{object}, or
#' if \code{object} is a \code{data.frame}, its row count.
#'
#' Two substitions are made \code{sprintf}-style. "\code{\%n}" is replaced with
#' the number of \code{object}, and "\code{\%c}" is replaced with the
#' comma-concatenated values of \code{object}, as in \code{\link{cc}}.
#'
#' \code{cn_and} uses \code{\link{cc_and}} instead of \code{cc}; \code{cn_or}
#' uses \code{\link{cc_or}}.
#'
#' @param object An n-vector, or \code{data.frame} with n rows.
#' @param singular The string to return if n = 1.
#' @param plural The string to return if n is in 0, 2, 3, 4, ...
#' @seealso \link{cc}
#' @name cn
NULL

cn_ <- function(FUN, n, object, singular, plural, ...) {
  msg <- ifelse(n == 1, singular, plural)
  msg <- gsub("%c", FUN(object), msg, fixed = TRUE)
  msg <- gsub("%n", n, msg, fixed = TRUE)
  msg
}

#' @rdname cn
#' @export
cn <- function(object, singular, plural = singular) {
  cn_(cc, length(object), object, singular, plural)
}

#' @rdname cn
#' @export
cn_and <- function(object, singular, plural = singular) {
  cn_(cc_and, length(object), object, singular, plural)
}

#' @rdname cn
#' @export
cn_or <- function(object, singular, plural = singular) {
  cn_(cc_or, length(object), object, singular, plural)
}

setGeneric("cn", signature = "object",
           function(object, singular, plural = singular) standardGeneric("cn"))
setGeneric("cn_and", signature = "object",
           function(object, singular, plural = singular) standardGeneric("cn_and"))
setGeneric("cn_or", signature = "object",
           function(object, singular, plural = singular) standardGeneric("cn_or"))

#' @rdname cn
#' @inheritParams cn
#' @import methods
#' @export
setMethod("cn", "data.frame",
          function(object, singular, plural = singular) {
            cn_(cc, nrow(object), object, singular, plural)
})

#' @rdname cn
#' @inheritParams cn
#' @import methods
#' @export
setMethod("cn_and", "data.frame",
          function(object, singular, plural = singular) {
            cn_(cc_and, nrow(object), object, singular, plural)
})

#' @rdname cn
#' @inheritParams cn
#' @import methods
#' @export
setMethod("cn_or", "data.frame",
          function(object, singular, plural = singular) {
            cn_(cc_or, nrow(object), object, singular, plural)
})
