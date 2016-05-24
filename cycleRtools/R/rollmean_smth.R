#' Rolling average smoothing.
#'
#' Smooth data with a right-aligned (zero-padded) rolling average.
#'
#' \code{rollmean_} is the core Rcpp function, which rolls over elements in
#' \code{x} by a window given in \code{window}; optionally applying exponential
#' weights and/or removing \code{NA}s. \code{rollmean_smth} is a wrapper for
#' \code{rollmean_} that only has a method for \code{cycleRdata} objects. The
#' latter will pre-process the data and permits what is effectively the
#' \code{window} argument being given in time units.
#'
#' @param data a dataset of class \code{cycleRdata}.
#' @param x numeric; values to be rolled over.
#' @param narm logical; should \code{NA}s be removed?
#' @param column the column name of the data to be smoothed, needn't be quoted.
#' @param smth.pd numeric; the time period over which to smooth (seconds).
#' @param window numeric; size of the rolling window in terms of elements in
#'   \code{x}.
#' @param ema logical; should the moving average be exponentially weighted?
#' @param deltat the sampling frequency of \code{data} in seconds per sample;
#'   typically 0.5 or 1. If \code{NULL}, this is estimated.
#' @param character.only are column name arguments given as character strings? A
#'   backdoor around non-standard evaluation.
#'
#' @return a vector of the same length as the \code{data[, column]}.
#'
#' @examples
#' \dontrun{
#' data(ridedata)
#'
#' ## Smooth power data with a 30 second moving average.
#' rollmean_smth(ridedata, power.W, 30)
#'
#' ## Or use an exponentially weighted moving average.
#' rollmean_smth(ridedata, power.W, 30, ema = TRUE)
#' }
#' @export
rollmean_smth <- function(data, column, smth.pd, deltat = NULL,
                          ema = FALSE, character.only = FALSE)
  UseMethod("rollmean_smth", data)
#' @export
rollmean_smth.default <- function(data, column, smth.pd, deltat = NULL,
                                  ema = FALSE, character.only = FALSE)
  format_error()
#' @export
rollmean_smth.cycleRdata <- function(data, column, smth.pd, deltat = NULL,
                                     ema = FALSE, character.only = FALSE) {
  if (!character.only)
    column <- as.character(substitute(column))

  data <- expand_stops(data, deltat, "timer.s")

  if (!is.null(attr(data, "new")))
    data[attr(data, "new"), "power.W"] <- 0 # Fill stops with zeros.

  out  <- rollmean_(data[, column], ema = ema,
                    window = (smth.pd / attr(data, "deltat")),
                    narm = TRUE)

  out[attr(data, "wo_expand")]
}
