#' GoldenCheetah (>v3.3) interface.
#'
#' Functions for interfacing R with
#' \href{https://github.com/GoldenCheetah/GoldenCheetah}{GoldenCheetah}.
#' Requires the \code{RCurl} package to be installed.
#'
#' As of GoldenCheetah (GC) version 3.3, the application is ran with a
#' background restful web service api to ease integration with external analysis software
#' (such as R). When an instance of GoldenCheetah is running, or the application
#' is initiated from the command line with the '--server' option, these
#' functions can be used to interface with athlete data. Relevant documentation
#' can be found
#' \href{https://github.com/GoldenCheetah/GoldenCheetah/wiki/UG_Special-Topics_REST-API-documentation}{here}.
#'
#' \code{GC_activity} behaves similarly to \code{\link{read_ride}} functions in
#' this package, importing data from saved GC .json files.
#'
#' \code{GC_metrics} returns summary metrics for either: all available rides if
#' \code{date.rng = NULL}; or rides within a specified date range if dates are
#' given.
#'
#' \code{GC_mmvs} retuns best maximal mean values for data specified in the
#' \code{type} argument. Possible options for \code{type} are: "watts", "hr",
#' "cad", "speed", "nm", "vam", "xPower", or "NP". See also \code{\link{mmv}}.
#'
#' @param athlete.name character; athlete of interest in the GoldenCheetah data
#'   directory. Typically of the form "First Last".
#' @param activity character; file path to a GoldenCheetah activity(.json) file.
#'   Typically located in "~/.goldencheetah/Athlete Name/activities/".
#' @param port http server port number. 12021 unless deliberatley changed in the
#'   httpserver.ini file.
#' @param format format activity data to an object of class "cycleRdata".
#'   Ensures compatibility with other functions in this package -- see
#'   \code{\link{read_ride}}.
#' @param date.rng a vector of length two that can be converted to an object of
#'   class \code{"Date"} via \code{\link{as.Date}}. Must be specified for
#'   \code{GC_mmvs}; optional for \code{GC_metrics}.
#' @param type the type of maximal mean values to return. See details.
#'
#' @name GC
NULL
# ------------------------------------------------------------------------------
GC_server_check <- function(URL) {
  msg <- paste("Localhost not found.",
               "Make sure GoldenCheetah is running and the listener is enabled.",
               sep = "\n")
  if (!RCurl::url.exists(URL))
    stop(msg, call. = FALSE)
}
#' @rdname GC
#' @export
GC_activity <- function(athlete.name, activity, port = 12021, format = TRUE) {
  URL          <- paste0("http://localhost:", port, "/") ; GC_server_check(URL)
  athlete.name <- gsub(" ", "%20", athlete.name)
  to_read <- paste0(
    URL,
    athlete.name,
    "/activity/",
    basename(activity),
    "?format=csv"
  )
  out <- read.csv(text = RCurl::getURL(to_read))
  if (format) {
    out        <- format_GC(out, filename = activity)
    class(out) <-  c("cycleRdata", "data.frame")
  }
  out
}
#' @rdname GC
#' @export
GC_metrics <- function(athlete.name, date.rng = NULL, port = 12021) {
  URL          <- paste0("http://localhost:", port, "/") ; GC_server_check(URL)
  athlete.name <- gsub(" ", "%20", athlete.name)
  to_read      <- paste0(URL, athlete.name)
  if (!is.null(date.rng)) {
    if (inherits(try(as.Date(date.rng), silent = TRUE), "try-error"))
      stop("Invalid date range.", call. = FALSE)
    date.rng <- format.Date(sort(date.rng), format = "%Y/%m/%d")
    to_read  <- paste0(to_read,
                       "?since=",  date.rng[1],
                       "&before=", date.rng[2])
  }
  out <- read.csv(text = RCurl::getURL(to_read))
  out
}
#' @rdname GC
#' @export
GC_mmvs <- function(type = "watts", date.rng = NULL, port = 12021) {
  if (is.null(date.rng))
    stop("date.rng argument is required.", call. = FALSE)
  if (inherits(try(as.Date(date.rng), silent = TRUE), "try-error"))
    stop("Invalid date range.", call. = FALSE)
  URL      <- paste0("http://localhost:", port, "/") ; GC_server_check(URL)
  # ----------------------------------------------------------------------
  date.rng <- format.Date(sort(date.rng), format = "%Y/%m/%d")
  to_read  <- paste0(URL,
                     "Jordan%20Mackie/meanmax/bests?",
                     "series=",  type,
                     "&since=",  date.rng[1],
                     "&before=", date.rng[2])
  out      <- read.csv(text = RCurl::getURL(to_read))
  out
}
