#' Get all global objects of an expression
#'
#' @param expr An R expression.
#' @param envir The environment from where to search for globals.
#' @param \dots Not used.
#' @param method A character string specifying what type of search algorithm to use.
#' @param tweak An optional function that takes an expression
#'        and returns a tweaked expression.
## @param dotdotdot A @character string specifying how to handle a
##        \emph{global} \code{\dots} if one is discovered.
#' @param substitute If TRUE, the expression is \code{substitute()}:ed,
#'        otherwise not.
#' @param mustExist If TRUE, an error is thrown if the object of the
#'        identified global cannot be located.  Otherwise, the global
#'        is not returned.
#' @param unlist If TRUE, a list of unique objects is returned.
#'        If FALSE, a list of \code{length(expr)} sublists.
#'
#' @return A \link{Globals} object.
#'
#' @details
#' There currently three methods for identifying global objects.
#'
#' The \code{"ordered"} search method identifies globals such that
#' a global variable preceeding a local variable with the same name
#' is not dropped (which the \code{"conservative"} method would).
#'
#' The \code{"conservative"} search method tries to keep the number
#' of false positive to a minimum, i.e. the identified objects are
#' most likely true global objects.  At the same time, there is
#' a risk that some true globals are not identified (see example).
#' This search method returns the exact same result as the
#' \code{\link[codetools]{findGlobals}()} function of the
#' \pkg{codetools} package.
#'
#' The \code{"liberal"} search method tries to keep the
#' true-positive ratio as high as possible, i.e. the true globals
#' are most likely among the identified ones.  At the same time,
#' there is a risk that some false positives are also identified.
#'
#' @example incl/globalsOf.R
#'
#' @seealso
#' Internally, the \pkg{\link{codetools}} package is utilized for
#' code inspections.
#'
#' @aliases findGlobals
#' @export
globalsOf <- function(expr, envir=parent.frame(), ..., method=c("ordered", "conservative", "liberal"), tweak=NULL, substitute=FALSE, mustExist=TRUE, unlist=TRUE) {
  method <- match.arg(method)

  if (substitute) expr <- substitute(expr)

  names <- findGlobals(expr, envir=envir, ..., method=method, tweak=tweak, substitute=FALSE, unlist=unlist)

  n <- length(names)
  needsDotdotdot <- (identical(names[n], "..."))
  if (needsDotdotdot) names <- names[-n]

  globals <- structure(list(), class=c("Globals", "list"))
  where <- list()
  for (name in names) {
    env <- where(name, envir=envir, inherits=TRUE)
    if (!is.null(env)) {
      where[[name]] <- env
      value <- get(name, envir=env, inherits=FALSE)
      if (is.null(value)) {
        globals[name] <- list(NULL)
      } else {
        globals[[name]] <- value
      }
    } else {
      where[name] <- list(NULL)
      if (mustExist) {
        stop(sprintf("Identified a global object via static code inspection (%s), but failed to locate the corresponding object in the relevant environments: %s", hexpr(expr), sQuote(name)))
      }
    }
  }

  if (needsDotdotdot) {
    if (exists("...", envir=envir, inherits=TRUE)) {
      where[["..."]] <- where("...", envir=envir, inherits=TRUE)
      ddd <- evalq(list(...), envir=envir, enclos=envir)
    } else {
      where["..."] <- list(NULL)
      ddd <- NA
    }
    class(ddd) <- c("DotDotDotList", class(ddd))
    globals[["..."]] <- ddd
  }

  attr(globals, "where") <- where

  globals
}
