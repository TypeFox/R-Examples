#' Smoothed data plot.
#'
#' Create a plot with both raw and smoothed data lines.
#'
#' @param data the dataset to be used.
#' @param x column identifier for the x axis data.
#' @param yraw column identifier for the (underlying) raw data.
#' @param ysmth column identifier for the smoothed data.
#' @param colour level identifier in \code{data} by which to colour lines. Or a
#'   colour name.
#' @param ... further arguments to be passed to \code{plot()}.
#' @param character.only are column name arguments given as character strings? A
#'   backdoor around non-standard evaluation.
#'
#' @examples
#' data(ridedata)
#'
#' ## Plot with a single blue line (default arguments):
#' smth_plot(ridedata, colour = "blue", main = "Single Colour",
#'           xlab = "Time (seconds)", ylab = "Power (watts)")
#'
#' ## Create some laps.
#' ridedata$lap <- ceiling(seq(from = 1.1, to = 5, length.out = nrow(ridedata)))
#' ## Plot with lap colours.
#' smth_plot(ridedata, timer.min, power.W, power.smooth.W, colour = "lap",
#'           xlab = "Time (mins)", ylab = "Power (watts)", main = "Lap Colours")
#'
#' @export
smth_plot <- function(data, x = "timer.s", yraw = "power.W", ysmth = "power.smooth.W",
                      colour = "lap", ..., character.only = FALSE) {
  if (!character.only) {
    x     <- as.character(substitute(x))
    yraw  <- as.character(substitute(yraw))
    ysmth <- as.character(substitute(ysmth))
  }
  # Fundamental plot.
  plot.default(x = data[, x], y = data[, yraw], type = "l", col = "gray", ...)
  # Overlay smooth line; handling colour argument.
  lwd <- 2
  if (!is.null(colour)) {
    # Argument is a colour name.
    if (!inherits(try(col2rgb(colour), silent = TRUE), "try-error")) {
      lines(x = data[, x], y = data[, ysmth], col = colour, lwd = lwd)
    } else if (!is.null(data[[colour]])) {
      # Argument is presumably a column name; does it exist?
      s <- unique(data[[colour]])
      for (i in s)
        lines(x = data[data[, colour] == i, x],
              y = data[data[, colour] == i, ysmth],
              col = base_pal(i, s),
              lwd = lwd)
    } else {
      # Column doesn't exist.
      lines(x = data[, x], y = data[, ysmth], col = colour, lwd = lwd)
    }
  } else {
    # NULL argument.
    lines(x = data[, x], y = data[, ysmth], col = "black", lwd = lwd)
  }
  invisible()
}
