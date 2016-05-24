#' Maximal mean values.
#'
#' Calculate maximal mean values for specified time periods.
#'
#' @param data a \strong{formatted} dataset produced by \code{read*()}.
#' @param column column in \code{data} giving the values of interest. Needn't be
#'   quoted.
#' @param windows window size(s) for which to generate best averages, given in
#'   seconds.
#' @param deltat the sampling frequency of \code{data} in seconds per sample;
#'   typically 0.5 or 1. If \code{NULL}, this is estimated.
#' @param character.only are column name arguments given as character strings? A
#'   backdoor around non-standard evaluation. Mainly for internal use.
#'
#' @return a matrix object with two rows: 1) best mean values and 2) the time
#'   at which those values were recorded
#'
#' @examples
#' data(ridedata)
#'
#' ## Best power for 5 and 20 minutes.
#' tsec <- c(5, 20) * 60
#' mmv(ridedata, power.W, tsec)
#'
#' ## Generate a simple critical power estimate.
#' tsec <- 2:20 * 60
#' pwrs <- mmv(ridedata, power.W, tsec)
#' m <- lm(pwrs[1, ] ~ {1 / tsec})  # Simple inverse model.
#' coef(m)[1]  # Intercept = critical power.
#'
#' ## More complex models...
#' m <- Pt_model(pwrs[1, ], tsec)
#' print(m)
#' ## Extract the asymptote of the exponential model.
#' coef(m)$exp["CP"]
#'
#' @seealso For a more generic and efficient version of this function, see
#'   \code{\link{mmv2}}
#'
#' @export
mmv <- function(data, column, windows, deltat = NULL,
                character.only = FALSE)
  UseMethod("mmv", data)
#' @export
mmv.default <- function(data, column, windows, deltat = NULL,
                        character.only = FALSE)
  format_error()
#'@export
mmv.cycleRdata <- function(data, column, windows, deltat = NULL,
                           character.only = FALSE) {

  if (!character.only)
    column <- as.character(substitute(column))

  data <- expand_stops(data, deltat, "timer.s")

  if (!is.null(attr(data, "new")))
    data[attr(data, "new"), "power.W"] <- 0 # Fill stops with zeros.

  mm_fn <- function(x, d, c) {
    mvavg <- rollmean_(d[, c], window = (x / attr(d, "deltat")),
                       ema = FALSE, narm = TRUE)
    best  <- max(mvavg, na.rm = TRUE)
    t     <- d[match(best, mvavg), "timer.s"]
    t     <- t - x # Start of the window rather than the end.
    c(best, t)
  }

  max_vals <- vapply(windows, FUN = mm_fn, d = data, c = column, numeric(2))
  rownames(max_vals) <- c("Best mean value", "Recorded @")
  colnames(max_vals) <- paste(windows)
  max_vals
}
