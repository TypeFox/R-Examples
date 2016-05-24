#' Wrapper for system2 with better return type and errorhandling.
#'
#' Wrapper for \code{\link{system2}} with better return type and errorhandling.
#'
#' @param command See \code{\link{system2}}.
#' @param args See \code{\link{system2}}.
#' @param stdout See \code{\link{system2}}.
#' @param stderr See \code{\link{system2}}.
#' @param wait See \code{\link{system2}}.
#' @param ... Further arguments passed to \code{\link{system2}}.
#' @param stop.on.exit.code [\code{logical(1)}]\cr
#'   Should an exception be thrown if an exit code greater 0 is generated?
#'   Can only be used if  \code{wait} is \code{TRUE}.
#'   Default is \code{wait}.
#' @return [\code{list}].
#'   \item{exit.code [integer(1)]}{Exit code of command. Given if wait is \code{TRUE}, otherwise \code{NA}. 0L means success. 127L means command was not found}
#'   \item{output [character]}{Output of command on streams. Only given is \code{stdout} or \code{stderr} was set to \code{TRUE}, otherwise \code{NA}.}
#' @export
system3 = function(command, args = character(0L), stdout = "", stderr = "", wait = TRUE, ..., stop.on.exit.code = wait) {
  if (stop.on.exit.code && !wait)
    stopf("stop.on.exit.code is TRUE but wait is FALSE!")
  output = NA_character_
  if (isTRUE(stdout) || isTRUE(stderr)) {
    wait = TRUE
    # here we wait anyway and output of cmd is returned
    ec = 0L
    suppressWarnings({
        withCallingHandlers({
            op = system2(command = command, args = args, stdout = stdout, stderr = stderr, wait = wait, ...)
          }, warning = function(w) {
            g = gregexpr("\\d+", w$message)[[1L]]
            start = tail(g, 1L)
            len = tail(attr(g, "match.length"), 1L)
            ec <<- as.integer(substr(w$message, start, start + len - 1L))
          })
      })
  } else {
    ec = system2(command = command, args = args, stdout = stdout, stderr = stderr, wait = wait, ...)
  }
  if (wait) {
    if (isTRUE(stdout) || isTRUE(stderr))
      output = op
  }
  if (stop.on.exit.code && ec > 0L) {
    args = collapse(args, " ")
    if (length(output) == 0L)
      output = ""
    else
      output = collapse(output, "\n")
    stopf("Command: %s %s; exit code: %i; output: %s", command, args, ec, output)
  }
  list(exit.code = ec, output = output)
}
