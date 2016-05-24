## time_utils.R
##   - Utility functions for dates and times
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Format time and data values into human-readable character vectors
##'
##' These functions convert date and time values into human-readable character
##' vectors.
##' \code{formatSeconds} formats time values given as a numerical vector denoting
##' seconds into human-readable character vectors, i.e. \code{formatSeconds(70)}
##' results in the string \code{"1 minute, 10 seconds"}.
##'
##' @param seconds A numeric vector denoting seconds.
##' @param secondDecimals The number of decimal places to show for seconds.
##'   Defaults to 2.
##' @return A character vector containg a human-readable representation of the
##'   given date/time.
##' @export
formatSeconds <- function(seconds, secondDecimals = 2) {
  formSecs <- function(seconds) {
    days <- floor(seconds / 86400)
    hours <- floor((seconds / 3600) - (days * 24))
    minutes <- floor((seconds / 60) - (days * 1440) - (hours * 60))
    restSeconds <- round(seconds %% 60, secondDecimals)
    # Coding Challenge:
    # If you know how to do exactly the following with sprintf, feel free to change it... ;-)
    paste(if (days == 1) paste(days, " day, ", sep = "") else if (days > 1) paste(days, " days, ", sep = "") else "",
          if (hours == 1) paste(hours, " hour, ", sep = "") else if (hours > 1) paste(hours, " hours, ", sep = "") else "",
          if (minutes == 1) paste(minutes, " minute, ", sep = "") else if (minutes > 1) paste(minutes, " minutes, ", sep = "") else "",
          if (restSeconds == 1) paste(restSeconds, " second", sep = "") else paste(restSeconds, " seconds", sep = ""),
          sep = "")
  }
  sapply(seconds, formSecs)
}
