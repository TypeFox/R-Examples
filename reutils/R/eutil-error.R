#' @include utils.R
#' @importFrom methods show
NULL

#' Class \code{"eutil_error"}:
#'
#' A container for handling errors when trying to parse XML files returned
#' by Entrez.
#' 
#' @field Error messages returned by the Entrez server when no search could be
#' performed like, e.g., "Invalid db name specified".
#' @field errmsg Error messages pertaining to the search like e.g.,
#' "PhraseNotFound".
#' @field wrnmsg Warnings like, e.g., "No items found."
#' 
#' @section Extends: All reference classes extend and inherit methods from
#'     \code{"\linkS4class{envRefClass}"}.
#' @seealso \code{\link{getError}}, \code{\linkS4class{eutil}}.
#' 
#' @name eutil_error-class
#' @keywords classes internal
#' @export
#' @examples
#' showClass("eutil_error")
eutil_error <- setRefClass(
  Class   = "eutil_error",
  fields  = c("error", "errmsg", "wrnmsg"),
  methods = list(
    initialize = function() {
      .self$error <- NULL
      .self$errmsg <- NULL
      .self$wrnmsg <- NULL
    },
    all_empty = function() {
      'Are all error fields \\code{NULL}?'
      is.null(error) && is.null(errmsg) && is.null(wrnmsg)
    },
    check_errors = function(object, verbose = TRUE) {
      'check if a \\code{linkS4class{eutil}} object contains errors'
      assertthat::assert_that(is(object, "eutil"))
      if (is.null(object$retmode()) || object$retmode() != "json") {
        x <- object$get_content("xml")
        .self$error <- xvalue(x, '//ERROR', default = NULL)
        if (verbose && !is.null(error)) {
          message('Error:\n\t', error)
        }
        errmsg_name  <- xname(x, '//ErrorList/*', default = NULL)
        .self$errmsg <- setNames(xvalue(x, '//ErrorList/*', default = NULL), errmsg_name)
        if (verbose && !is.null(errmsg)) {
          message('Error(s):\n\t', 
                  paste(paste(names(errmsg), errmsg, sep = "\t"), collapse = "\n\t"))
        }
        wrnmsg_name  <- xname(x, '//WarningList/*', default = NULL)
        .self$wrnmsg <- setNames(xvalue(x, '//WarningList/*', default = NULL), wrnmsg_name) 
        if (verbose && !is.null(wrnmsg)) {
          message('Warning(s):\n\t', 
                  paste(paste(names(wrnmsg), wrnmsg, sep = "\t"), collapse = "\n\t"))
        }
      }
    },
    show = function() {
      if (all_empty()) {
        cat("No errors", sep = "\n")
        invisible(.self)
      } else {
        error  %&&% methods::show(error)
        errmsg %&&% methods::show(errmsg)
        wrnmsg %&&% methods::show(wrnmsg)
        invisible(.self)
      }
    }
  )
)
