#' Plot gap-filled dendrometer series
#'
#' @description The function creates a plot with gap-filled and original dendrometer series.
#'
#' @usage fill_plot(dm.data, dm.gpf, sensor = 1, year = NULL,
#'           period = NULL)
#'
#' @param dm.data a \code{data.frame} with a timestamp (\code{\%Y-\%m-\%d \%H:\%M:\%S} format) as row names, and dendrometer series in columns. Output as created using code from the \code{Import dendrometer data} vignette.
#' @param dm.gpf a \code{data.frame} with gap-filled dendrometer series as produced by \code{\link{fill_gaps}}.
#' @param sensor a \code{numeric} specifying the sensor to be plotted (by column number). Defaults to 1 (first dendrometer series in both \code{data.frames}).
#' @param year a \code{numeric} specifying the year(s) to be plotted. Defaults to the first year in the dataset. Two consecutive years (e.g., for a growing season at the Southern Hemisphere) can be defined with \code{year = c(year1, year2)}.
#' @param period a \code{numeric} indicating the period to be plotted, specified using day of year values (begin and end). Defaults to the complete data period.
#'
#' @details The function creates a plot showing the gap-filling results for a single dendrometer series over a specified time window. Although the function is intended to plot short time periods (within a growing season), it can plot two calendar years at maximum (e.g., 2014-2015), thereby allowing the visualization of a complete growing season at the Southern Hemisphere as well.
#'
#' @return Plot.
#'
#' @author Olivier Bouriaud, Ernst van der Maaten and Marieke van der Maaten-Theunissen.
#'
#' @examples
#' \dontrun{
#'
#' data(dmCD)
#' ## creating some artificial gaps (for demonstration purposes):
#' dmCD[c(873:877,985:990),1] <- NA
#' dm.gpf <- fill_gaps(dmCD, Hz = 0.01)
#' fill_plot(dmCD, dm.gpf, period = c(137,144))
#' }
#'
#' @import zoo
#' @import pspline
#' @import graphics
#'
#' @export fill_plot
#'
fill_plot <- function(dm.data, dm.gpf, sensor = 1, year = NULL, period = NULL)
{
  nm1 <- deparse(substitute(dm.data))
  nm2 <- deparse(substitute(dm.gpf))

  if(!is.dendro(dm.data) | !is.dendro(dm.gpf)) {
    stop("there is a problem with the input data")
  }
  if(nrow(dm.data) != nrow(dm.gpf)) {
    stop(paste("'", nm1, "' and '", nm2,"' have different numbers of rows", sep = ""))
  }
  if(ncol(dm.data) != ncol(dm.gpf)) {
    stop(paste("'", nm1, "' and '", nm2,"' have different numbers of columns", sep = ""))
  }
  if(ncol(dm.gpf) == 1 && sensor > ncol(dm.gpf)) {
    sensor <- 1
    warning(paste("'", nm1, "' contains only one dendrometer series, 'sensor' was set to 1", sep = ""))
  }
  if(ncol(dm.gpf) > 1 && sensor > ncol(dm.gpf)) {
    stop(paste("'sensor' should be between", 1, "and", ncol(dm.gpf), sep = " "))
  }

  resolution1 <- dendro.resolution(dm.data, "hours")
  resolution2 <- dendro.resolution(dm.gpf, "hours")
  if(resolution1 != resolution2) {
    stop(paste("'", nm1, "' and '", nm2, "' have a different temporal resolution", sep = ""))
  }
  period.var <- 24 / resolution1

  dm.data$YEAR <- as.numeric(strftime(rownames(dm.data), format = "%Y"))
  dm.gpf$YEAR <- as.numeric(strftime(rownames(dm.gpf), format = "%Y"))

  dm.data$DOY <- as.numeric(strftime(rownames(dm.data), format = "%j"))
  dm.gpf$DOY <- as.numeric(strftime(rownames(dm.gpf), format = "%j"))

  min.yr <- min(dm.gpf$YEAR)
  max.yr <- max(dm.gpf$YEAR)

  if(is.null(year)) {
    year <- min.yr
  }
  if(length(year) == 2 && year[1] == year[2]) {
    year <- year[1]
  }
  if(length(year) == 1 && min.yr == max.yr && (year < min.yr || year > max.yr)) {
    stop(paste("'year' should be", min.yr, sep = " "))
  }
  if(length(year) == 1 && min.yr != max.yr && (year < min.yr || year > max.yr)) {
    stop(paste("'year' should be between", min.yr, "and", max.yr, sep = " "))
  }
  if(length(year) == 2 && year[1] > year[2]) {
    stop("first entry for 'year' should be lower than second entry")
  }
  if(length(year) == 2 && year[2]-1 != year[1]) {
    stop("two consecutive years (or a single year) should be selected for 'year'")
  }
  if(length(year) == 2 && min.yr == max.yr && (TRUE %in% (year < min.yr) || TRUE %in% (year > max.yr))) {
    stop(paste("'year' should be", min.yr, sep = " "))
  }
  if(length(year) == 2 && (TRUE %in% (year < min.yr) || TRUE %in% (year > max.yr))) {
    stop(paste("'year' should be between", min.yr, "and", max.yr, sep = " "))
  }

  if(length(year) == 1) {
    sub.data <- dm.data[dm.data$YEAR == year,]
    sub.gpf <- dm.gpf[dm.gpf$YEAR == year,]
    min.doy <- min(sub.gpf$DOY)
    max.doy <- max(sub.gpf$DOY)
  }

  if(length(year) == 2) {
    sub.data <- dm.data[dm.data$YEAR == year[1] | year[2],]
    sub.gpf <- dm.gpf[dm.gpf$YEAR == year[1] | year[2],]
    min.doy <- min(sub.gpf[sub.gpf$YEAR == year[1], "DOY"])
    max.doy <- max(sub.gpf[sub.gpf$YEAR == year[2], "DOY"])
  }

  if(!is.null(period) && length(period) != 2) {
    warning("'period' should have two members; the default (NULL) is used")
    period <- NULL
  }
  if(!is.null(period) && !is.numeric(period)) {
    warning("'period' should identify a period based on numeric DOY values")
    period <- NULL
  }
  if(is.null(period)) {
    doy.start <- min.doy
    doy.end <- max.doy
  } else {
      doy.start <- min(period)
      doy.end <- max(period)
  }

  if(length(year) == 1 && (doy.start < min.doy || doy.end > max.doy)) {
    stop(paste("DOY values in 'period' should be between", min.doy, "and", max.doy, sep = " "))
  }
  if(length(year) == 2 && (doy.start < min.doy || doy.end > max.doy)) {
    stop(paste("DOY values in 'period' should be between", min.doy, "(year 1) and", max.doy, "(year 2)", sep = " "))
  }

  if(length(year) == 1) {
    sub.data2 <- sub.data[sub.data$DOY >= doy.start & sub.data$DOY <= doy.end, ]
    sub.gpf2 <- sub.gpf[sub.gpf$DOY >= doy.start & sub.gpf$DOY <= doy.end, ]
  }
  if(length(year) == 2) {
    sub.data2 <- sub.data[(sub.data$YEAR == year[1] & sub.data$DOY >= doy.start) |
                           (sub.data$YEAR == year[2] & sub.data$DOY <= doy.end), ]
    sub.gpf2 <- sub.gpf[(sub.gpf$YEAR == year[1] & sub.gpf$DOY >= doy.start) |
                           (sub.gpf$YEAR == year[2] & sub.gpf$DOY <= doy.end), ]
  }

  ts.data <- ts(zoo(sub.data2[, sensor], rownames(sub.data2)), start = doy.start, frequency = period.var)
  ts.gpf <- ts(zoo(sub.gpf2[, sensor], rownames(sub.gpf2)), start = doy.start, frequency = period.var)

  if(length(year) == 1) {
    plot(ts.gpf, type = "l", col = "darkorange1", xlab = "Day of year", ylab = "Stem-size variation",
         las = 1, main = paste("Sensor ", colnames(sub.data2)[sensor], " (", year, ")",
                               sep = ""), xlim = range(pretty(c(doy.start, doy.end))))
    lines(ts.data, lwd = 1.5)
  }

  if(length(year) == 2) {
    if((year[1]%%4 == 0) & ((year[1]%%100 != 0) | (year[1]%%400 == 0))){
    plot(ts.gpf, type = "l", col = "darkorange1", xlab = "Day of year", ylab = "Stem-size variation",
         las = 1, main = paste("Sensor ", colnames(sub.data2)[sensor], " (", year[1], "-", year[2], ")",
                              sep = ""), xaxt = "n")
    doy <- pretty(c(doy.start:366, (366+1):(366+doy.end)))
    doy.lab <- doy%%367+1
    doy.lab[which(doy < 367)] <- doy.lab[which(doy < 367)]-1
    axis(1, at = doy, labels = doy.lab)
    lines(ts.data, lwd = 1.5)
    }
    else {
      plot(ts.gpf, type = "l", col = "darkorange1", xlab = "Day of year", ylab = "Stem-size variation",
           las = 1, main = paste("Sensor ", colnames(sub.data2)[sensor], " (", year[1], "-", year[2], ")",
                                 sep = ""), xaxt = "n")
      doy <- pretty(c(doy.start:365, (365+1):(365+doy.end)))
      doy.lab <- doy%%366+1
      doy.lab[which(doy < 366)] <- doy.lab[which(doy < 366)]-1
      axis(1, at = doy, labels = doy.lab)
      lines(ts.data, lwd = 1.5)
    }
  }
}
