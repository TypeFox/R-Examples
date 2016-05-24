# luzlogr - user-visible logging functions

# -----------------------------------------------------------------------------
#' Open a new logfile
#'
#' @param file Name of logfile (character or writeable \code{\link{connection}})
#' @param loglevel Minimum priority level (numeric, optional)
#' @param append Append to logfile? (logical, optional)
#' @param sink Send all console output to logfile? (logical, optional)
#' @return Invisible fully-qualified name of log file.
#' @details Open a new logfile. Messages will only appear in the logfile
#' if their \code{level} exceeds the log's \code{loglevel};
#' this allows you to easily change the amount of detail being logged.
#'
#' Re-opening a logfile will erase the previous output unless \code{append}
#' is TRUE. Opening a new logfile when one is already open will temporarily
#' switch logging to that new file.
#'
#' If \code{sink} is TRUE, all screen output will be captured (via \code{\link{sink}}).
#' @examples
#' logfile <- openlog("test.log")
#' printlog("message")
#' closelog()
#' readLines(logfile)
#' @export
#' @seealso \code{\link{printlog}} \code{\link{closelog}}
openlog <- function(file, loglevel = -Inf, append = FALSE, sink = FALSE) {

  # Sanity checks
  assert_that(is.numeric(loglevel))
  assert_that(is.logical(append))
  assert_that(is.logical(sink))

  if(is.character(file)) {
    file <- file(file)  # THAT'S confusing! Change text filename to connection
  }

  if(inherits(file, "connection")) {  # connection
    closeit <- !isOpen(file)
    if(!isOpen(file)) {
      msg("Opening connection")
      open(file, if(append) "a" else "w")
    }
    description <- summary(file)$description
  }
  else stop("'file' must be a character string or a connection")

  # Open a new sink, if requested
  if(sink) {
    sink(file, split = TRUE, append = append)
  }

  # Create a new log in our internal data structure
  newlog(logfile = file, loglevel = loglevel, sink = sink,
         description = description, closeit = closeit)

  printlog("Opening", description, level = Inf)
  invisible(normalizePath(description))
} # openlog

# -----------------------------------------------------------------------------
#' Log a message
#'
#' @param ... Expressions to be printed to the log
#' @param level Priority level (numeric, optional)
#' @param ts Print preceding timestamp? (logical, optional)
#' @param cr Print trailing newline? (logical, optional)
#' @param flag Flag this message (e.g. error or warning) (logical, optional)
#' @param flush Immediately flush output to file (logical, optional)
#' @return Invisible success (TRUE) or failure (FALSE).
#' @details Logs a message, which consists of zero or more printable objects.
#' Simple objects (numeric and character) are printed together on a single
#' line, whereas complex objects (data frames, etc) start on a new line by
#' themselves.
#'
#' If the current log was opened with \code{sink} = TRUE,
#' messages are printed to the screen, otherwise not. Messages can be flagged;
#' \code{flaglog} assumes
#' that the message is to be flagged, whereas \code{printlog} does not.
#'
#' Messages will only appear in the logfile if their \code{level} exceeds
#' the log's \code{loglevel}; this allows you to easily change the amount of
#' detail being logged.
#' @note A message's preceding timestamp and following carriage return can be
#' suppressed using the \code{ts} and \code{cr} parameters.
#' @examples
#' logfile <- openlog("test.log")
#' printlog("message")
#' printlog(1, "plus", 1, "equals", 1 + 1)
#' closelog()
#' readLines(logfile)
#'
#' logfile <- openlog("test", loglevel = 1)
#' printlog("This message will not appear", level = 0)
#' printlog("This message will appear", level = 1)
#' closelog()
#' readLines(logfile)
#' @export
#' @seealso \code{\link{openlog}} \code{\link{closelog}}
printlog <- function(..., level = 0, ts = TRUE, cr = TRUE, flag = FALSE, flush = FALSE) {

  # Sanity checks
  assert_that(is.numeric(level))
  assert_that(is.logical(ts))
  assert_that(is.logical(cr))
  assert_that(is.logical(flag))
  assert_that(is.logical(flush))

  args <- list(...)

  # Make sure there's an open log file available
  loginfo <- getloginfo()
  if(is.null(loginfo)) {
    return(invisible(FALSE))
  }

  msg("Writing to ", loginfo$description)

  # If someone has messaged with sink(), by opening or closing
  # sinks, things will get screwed up. Warn in this case.
  if(sink.number() != loginfo$sink.number) {
    warning("Current sink number doesn't match one in log data")
  }

  # Messages are only printed if their level exceeds the log's level (or an error)
  if(level >= loginfo$loglevel | flag) {
    if(loginfo$sink) { # If capturing everything, output to screen
      file <- stdout()
    } else {  # otherwise, file
      file <- loginfo$logfile
    }

    # Print a special message if warning (flag) condition
    if(flag) {
      setlogdata("flags", loginfo$flags + 1)
      flagmsg <- "** Flagged message: **\n"
      #       if(loginfo$sink) {
      #         message(flagmsg)
      #       }
      cat(flagmsg, file = file, append = TRUE)
    }

    # Print a timestamp...
    if(ts) cat(date(), " ", file = file, append = TRUE)

    # ...and then the object(s)
    for(i in seq_along(args)) {
      x <- args[[i]]
      # simple objects are printed together on a line
      if(mode(x) %in% c("numeric", "character")) {
        cat(x, " ", file = file, append = TRUE, sep = "")
      } else { # more complex; let print() handle it
        cat("\n", file = file, append = TRUE)
        if(loginfo$sink) {
          print(x)
        } else {
          utils::capture.output(x, file = file, append = TRUE)
        }
      }
    }

    if(cr) cat("\n", file = file, append = TRUE)

    if(flush) {
      flush(file)   # force output to logfile
    }
  } else {
    msg("Message level not high enough")
  }

  invisible(TRUE)
} # printlog

# -----------------------------------------------------------------------------
#' @rdname printlog
#' @export
flaglog <- function(...) printlog(..., flag = TRUE)

# -----------------------------------------------------------------------------
#' Close current logfile
#'
#' @param sessionInfo Append \code{\link{sessionInfo}} output? (logical, optional)
#' @return Number of flagged messages (numeric).
#' @details Close current logfile. The number of flagged messages is returned,
#' invisibly. Note that if \code{options(luzlogr.close_on_error = TRUE)} is set, then
#' if an error occurs, all log files will be automatically closed. This behavior
#' is not currently enabled by default.
#'
#' Logs are stored on a stack, and so when one is closed, logging
#' output returns to the previous log (if any).
#' @note If the log was being written to a \code{\link{connection}},
#' \code{closelog} will return the connection to its pre-logging state,
#' whether open or closed.
#' @examples
#' logfile <- openlog("A.log")
#' printlog("message to A", flag = TRUE)
#' logfile <- openlog("B.log")
#' printlog("message to B")
#' flagcountB <- closelog()
#' flagcountA <- closelog(sessionInfo = FALSE)
#' @export
#' @seealso \code{\link{openlog}} \code{\link{printlog}}
closelog <- function(sessionInfo = TRUE) {

  # Make sure there's an open log file available to close
  loginfo <- getloginfo()
  if(is.null(loginfo)) return(invisible(NULL))

    description <- summary(loginfo$logfile)$description

  printlog("Closing", description, "flags =", loginfo$flags, level = Inf)

  # Remove sink, if applicable
  if(loginfo$sink & sink.number()) sink()

  # Append sessionInfo() to file
  if(sessionInfo) {
    cat("-------\n", file = loginfo$logfile, append = TRUE)
    utils::capture.output(sessionInfo(), file = loginfo$logfile, append = TRUE)
  }

  # Close file or connection, if necessary
  if(loginfo$closeit) {
    close(loginfo$logfile)
  }

  # Remove log from internal data structure
  removelog()

  invisible(loginfo$flags)
} # closelog
