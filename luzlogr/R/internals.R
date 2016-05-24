# luzlogr - backend support functions
#
# Keep a clean separation between implementation of loginfo data and use

PKG.ENV <- new.env()    # environment in which to store logging info
LOGINFO <- ".loginfo"   # name of storage variable

# LOGINFO, above, is implemented as a list of lists.
# The first-order list operates as a stack, where the last entry holds the
# currently-active log information. Each entry, in turn, is a list of
# information about that particular log. `newlog` pushes a new entry on to
# the stack, while `closelog` pops the last one off.

# -----------------------------------------------------------------------------
#' Create new log
#'
#' @param logfile Name of log file (character or connection)
#' @param loglevel Minimum priority level (numeric)
#' @param sink Send all console output to logfile? (logical)
#' @param description Description (character)
#' @param closeit File should be closed when log closes? (logical)
#' @details This handles internal data tracking only, not the file on disk.
#' @keywords internal
newlog <- function(logfile, loglevel, sink, description, closeit) {

  # Sanity checks
  assert_that(inherits(logfile, "connection"))
  assert_that(is.numeric(loglevel))
  assert_that(is.logical(sink))
  assert_that(is.character(description))
  assert_that(is.logical(closeit))

  # If log data structure doesn't exist, create a (hidden) variable
  # in the package environment to store it
  if(!exists(LOGINFO, envir = PKG.ENV)) {
    msg("Creating empty log data structure")
    assign(LOGINFO, list(), envir = PKG.ENV)
  }

  # Create a new entry for this particular log-open request
  newloginfo <- list(logfile = logfile,
                     loglevel = loglevel,
                     sink = sink,
                     description = description,
                     closeit = closeit,
                     sink.number = sink.number(),
                     flags = 0)

  # Does this connection already exist?
  # TODO

  loginfo <- get(LOGINFO, envir = PKG.ENV)
  loginfo[[length(loginfo) + 1]] <- newloginfo
  assign(LOGINFO, loginfo, envir = PKG.ENV)
  msg("Logs: ", length(loginfo))
} # newlog

# -----------------------------------------------------------------------------
#' Remove current log
#'
#' @return The just-removed log info.
#' @details This handles internal data tracking only, not the file on disk.
#' @keywords internal
removelog <- function() {
  # Get the current log data
  if(exists(LOGINFO, envir = PKG.ENV)) {
    loginfo <- get(LOGINFO, envir = PKG.ENV)

    if(length(loginfo)) {
      loginfotop <- loginfo[[length(loginfo)]]
      loginfo[[length(loginfo)]] <- NULL
      assign(LOGINFO, loginfo, envir = PKG.ENV)
      msg("Logs: ", length(loginfo))
      return(loginfotop)
    }
  }

  warning("No log available")
  NULL
} # removelog

# -----------------------------------------------------------------------------
#' Get current log info
#'
#' @return A list with information about current (active) log.
#' @details This handles internal data tracking only, not the file on disk.
#' @keywords internal
getloginfo <- function() {
  # Get the current log data
  if(exists(LOGINFO, envir = PKG.ENV)) {
    loginfo <- get(LOGINFO, envir = PKG.ENV)

    if(length(loginfo)) {
      return(loginfo[[length(loginfo)]])
    }
  }
  warning("No log available")
  NULL
} # getloginfo

# -----------------------------------------------------------------------------
#' Set log data
#'
#' @param logdata Name of datum to set (character)
#' @param value Value
#' @details This handles internal data tracking only, not the file on disk.
#' @keywords internal
setlogdata <- function(datum, value) {

  assert_that(is.character(datum))

  # Get the current log data
  if(exists(LOGINFO, envir = PKG.ENV)) {
    loginfo <- get(LOGINFO, envir = PKG.ENV)

    if(length(loginfo)) {
      loginfotop <- loginfo[[length(loginfo)]]
      assert_that(datum %in% names(loginfotop))

      # Currently we only allow changing the 'flags' value
      if(datum == "flags") {
        loginfotop$flags <- value
      } else {
        stop("Error: can't modify", datum)
      }

      loginfo[[length(loginfo)]] <- loginfotop
      assign(LOGINFO, loginfo, envir = PKG.ENV)
      return(invisible(NULL))
    }
  }

  stop("No log available")
} # setlogdata

# -----------------------------------------------------------------------------
#' Return number of current logs
#'
#' @return Number of current logs (numeric).
#' @keywords internal
nlogs <- function() {
  if(exists(LOGINFO, envir = PKG.ENV)) {
    length(get(LOGINFO, envir = PKG.ENV))
  } else
    0
} # nlogs
