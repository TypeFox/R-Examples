pwr_tf <- function(x, e = 4) {mean(x ^ e, na.rm = TRUE) ^ (1 / e)}

#' Summary metrics.
#'
#' Common summary measures of interest to cyclists.
#'
#' \code{NP} calculates a Normalised Power value. "Normalised Power" is a
#' registered trademark of Peaksware Inc.
#'
#' \code{xPower}; Dr. Philip Skiba/Golden Cheetah's answer to NP.
#'
#' \code{pwr_TRIMP}: Power-Based TRaining IMPulse. Calculates a
#' \emph{normalised} TRIMP value using power data. This is a power-based
#' adaptation of Bannister's TRIMP, whereby critical power (CP) is assumed to
#' represent 90% of heart rate ratio (HRR). The final TRIMP score is normalised
#' to the score associated with one-hour's riding at CP, to aid interpretation.
#'
#' \code{ride_time} is a simple function for calculating ride time, as opposed
#' to elapsed time.
#'
#' \code{TSS} calculates a Training Stress Score (TSS). TSS is a registered
#' trademark of Peaksware Inc.
#'
#' @param data a "cycleRdata" object, produced from a \code{\link{read_ride}}
#'   function.
#' @param CP a Critical Power value - e.g. CP or FTP.
#' @param x a vector of time values.
#' @param deltat numeric; the typical interval between time values, if
#'   \code{NULL} a best estimate is used.
#'
#' @references Morton, R.H., Fitz-Clarke, J.R., Banister, E.W., 1990. Modeling
#'   human performance in running. Journal of Applied Physiology 69, 1171-1177.
#'
#' @return a single numeric value.
#' @name summary_metrics
#'
#' @examples
#' data(ridedata)
#'
#' ## Display all summary metrics with an *apply call.
#' fns   <- list("ride_time", "xPower", "NP", "pwr_TRIMP", "TSS")
#' argl  <- list(data = ridedata, x = ridedata$timer.s, CP = 300)
#' metrs <- vapply(fns, function(f) {
#'   do.call(f, argl[names(argl) %in% names(formals(f))])
#' }, numeric(1))
#'
#' names(metrs) <- fns
#' print(metrs)
# ------------------------------------------------------- #
#' @export
#' @rdname summary_metrics
ride_time <- function(x, deltat = NULL) {
  if (is.null(deltat)) {
    diffs <- table(Diff(x))
    diffs <- rbind(count = unname(diffs), val = as.numeric(names(diffs)))
    # deltat should represent > 75% of sample time differences.
    if ((diffs["count", which.max(diffs["count", ])] / length(x)) < 0.75)
      stop(bad_deltat_msg())
    deltat <- diffs["val", which.max(diffs["count", ])]
  }
  full <- seq(from = 0, to = max(x, na.rm = TRUE), by = deltat)
  dt   <- c(0, Diff(full))
  rt   <- sum(dt[full %in% x])
  rt
}
#' @export
#' @rdname summary_metrics
xPower <- function(data) {
  if (!is.cycleRdata(data)) format_error()
  pwr_tf(rollmean_smth(data, "power.W", 25, ema = TRUE), 4)
}
#' @export
#' @rdname summary_metrics
NP <- function(data) {
  if (!is.cycleRdata(data)) format_error()
  pwr_tf(rollmean_smth(data, "power.W", 30, ema = FALSE), 4)
}
#' @export
#' @rdname summary_metrics
pwr_TRIMP <- function(data, CP = attr(data, "CP")) {
  if (!is.cycleRdata(data)) format_error()
  pwr  <- pwr_tf(rollmean_smth(data, "power.W", 25, ema = TRUE), 4)
  tmin <- ride_time(data$timer.s) / 60
  int  <- {pwr / CP}
  ((tmin * 0.9 * int * exp(1.92 * 0.9 * int)) /
    (60 * 0.9 * exp(1.92 * 0.9))) * 100
}
#' @export
#' @rdname summary_metrics
TSS <- function(data, CP = attr(data, "CP")) {
  if (!is.cycleRdata(data)) format_error()
  timeTotal <- ride_time(data$timer.s)
  NormPwr <- NP(data)
  ((timeTotal * NormPwr * (NormPwr / CP)) /
    (CP * (60 ^ 2))) * 100
}
