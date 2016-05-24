#' Lightweight logging for R
#'
#' This package provides flexible but lightweight logging facilities for R scripts.
#' Supports priority levels for logs and messages, flagging messages,
#' capturing script output, switching logs, and logging to files or connections.
#'
#' @import assertthat
#' @docType package
#' @name luzlogr
NULL


DEBUG <- FALSE
msg <- function(...) if(DEBUG) message(...)


# On load, set up default package behavior
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.luzlogr <- list(
    luzlogr.close_on_error = FALSE
  )

  toset <- !(names(op.luzlogr) %in% names(op))
  if(any(toset)) {
    options(op.luzlogr[toset])
  }

  invisible()
} # .onLoad


# On attach, set up error handler
.onAttach <- function(libname, pkgname) {
  if(getOption("luzlogr.close_on_error")) {

    # Save whatever is currently assigned as an error handler (possibly NULL)
    assign(".errorfunc", getOption("error"), envir = PKG.ENV)

    closeall <- function() {
      msg("Closing all open logs")
      while(nlogs() > 0) {
        closelog(sessionInfo = FALSE)
      }
      # Execute whatever was handling error before us
      eval(get(".errorfunc", envir = PKG.ENV))
    }

    # Put in our new error handler
    options(error = closeall)
  }
} # .onAttach
