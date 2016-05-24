#' Split seconds into handy chunks of time.
#'
#' Note that a year is simply defined as exactly 365 days.
#'
#' @param seconds [\code{numeric(1)}]\cr
#'   Number of seconds. If not an integer, it is rounded down.
#' @param unit [\code{character(1)}]\cr
#'   Largest unit to split seconds into.
#'   Must be one of: \code{c("years", "days", "hours", "minutes", "seconds")}.
#'   Default is \dQuote{years}.
#' @return [\code{numeric(5)}]. A named vector containing the
#' \dQuote{years}, \dQuote{days}, \dQuote{hours}, \dQuote{minutes}
#' and \dQuote{seconds}. Units larger than the given \code{unit} are
#' \code{NA}.
#' @export
#' @examples
#' splitTime(1000)
splitTime = function(seconds, unit = "years") {
  assertNumber(seconds)
	assertChoice(unit, c("years", "days", "hours", "minutes", "seconds"))
  divider = c(31536000L, 86400L, 3600L, 60L, 1L)
  res = setNames(rep.int(NA_integer_, 5L),
                 c("years", "days", "hours", "minutes", "seconds"))
  start = which(names(res) == unit)
  for (i in start:length(divider)) {
    res[i] = seconds %/% divider[i]
    seconds = seconds - res[i] * divider[i]
  }
  ## Make sure all values are integral and do _not_ strip names:
  viapply(res, as.integer)
}
