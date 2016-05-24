#' Simple logger which outputs to a file.
#'
#' Creates a simple file logger closure to log to a file, including time stamps.
#' An optional buffer holds the last few log messages.
#'
#' @param logfile [\code{character(1)}]\cr
#'   File to log to.
#' @param touch [\code{logical(1)}]\cr
#'   Should the file be created before the first log message?
#'   Default is \code{FALSE}.
#' @param keep [\code{integer(1)}]\cr
#'   Number of log messages to keep in memory for quick access.
#'   Default is \code{10}.
#' @return [\code{\link{SimpleFileLogger}}]. A list with following functions:
#'   \item{log [\code{function(msg)}]}{Send log message.}
#'   \item{getMessages [\code{function(n)}]}{Get last \code{n} log messages.}
#'   \item{clear [\code{function()}]}{Resets logger and deletes log file.}
#'   \item{getSize [\code{function()}]}{Returns the number of logs written.}
#'   \item{getLogfile [\code{function()}]}{Returns the full file name logs are written to.}
#' @export
#' @aliases SimpleFileLogger
makeSimpleFileLogger = function(logfile, touch = FALSE, keep = 10L) {
  assertString(logfile)
  assertFlag(touch)
  keep = asCount(keep)
  assertDirectory(dirname(logfile), "w")
  if (touch && !file.create(logfile))
    stopf("Could not create file '%s'", logfile)
  if (keep)
    buffer = circularBuffer("character", keep)
  n.lines = 0L

  makeS3Obj("SimpleFileLogger",
    log = function(msg) {
      if (keep)
        buffer$push(msg)
      if (!touch && n.lines == 0L && !file.create(logfile))
        stopf("Could not create file '%s'", logfile)
      catf("<%s> %s", as.character(Sys.time()), msg, file = logfile, append = TRUE, newline = TRUE)
      n.lines <<- n.lines + 1L
    },
    getMessages = function(n) {
      if (!keep || n > keep)
        return(sub("^<.+> (.*)", "\\1", tail(readLines(logfile), n)))
      buffer$get(n)
    },
    clear = function() {
      if (keep)
        buffer$clear()
      n.lines <<- 0L
      file.remove(logfile)
    },
    getSize = function() {
      n.lines
    },
    getLogfile = function() {
      logfile
    }
  )
}

circularBuffer = function(type, capacity) {
  st = vector(type, capacity)
  stored = 0L
  pos = 0L

  list(
    push = function(x) {
      pos <<- pos %% capacity + 1L
      stored <<- min(capacity, stored + 1L)
      st[[pos]] <<- x
    },
    get = function(n = 1L) {
      head(st[rev((seq_len(stored) + pos - 1L) %% stored + 1L)], n)
    },
    stored = function() {
      stored
    },
    clear = function() {
      st <<- vector(type, capacity)
      stored <<- 0
      pos <<- 0
    }
  )
}
