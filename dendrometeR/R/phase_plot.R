#' Plot stem-cyclic phases
#'
#' @description The function creates a plot showing the three distinct phases of contraction, expansion and stem-radius increment (Deslauriers et al. 2011) for dendrometer series from a \code{data.frame} as produced by \code{\link{phase_def}}.
#'
#' @usage phase_plot(dm.gpf, dm.phase, sensor = NULL, period = NULL, colPhases = NULL, ...)
#'
#' @param dm.gpf a \code{data.frame} with gap-filled dendrometer series as produced by \code{\link{fill_gaps}}.
#' @param dm.phase a \code{data.frame} with numbers indicating the different stem-cyclic phases. Output of \code{\link{phase_def}}.
#' @param sensor a \code{numeric} specifying the sensor to be plotted (by column number). Alternatively, \code{sensor} can be a \code{character} with column names. Concatenations and sequences are allowed for plotting phase definitions of multiple sensors at once. Defaults to all sensors in \code{dm.gpf} and \code{dm.phase}.
#' @param period a \code{numeric} indicating the period to be plotted, specified using day of year values (begin and end). Defaults to the complete data period. Alternatively, \code{period} can be a \code{character} of two time stamps, indicating the begin and end date of the period to be plotted.
#' @param colPhases a \code{vector} of length 3, specifying custom colors to be used for the three stem-cyclic phases. Defaults to the first three colors from the current \code{\link{palette}}.
#' @param ... additional graphical parameters (see \code{\link{par}}).
#'
#' @return Plot showing stem-cyclic phases on dendrometer series.
#'
#' @details The function plots phases of contraction, expansion and stem-radius increment along (one or more) dendrometer series. If more series are plotted (default), colors for the different lines can be defined using the \code{col} argument for graphical devices (see \code{\link{par}}). Note: if there are not enough custom colors, the function will repeat the last one used. If no colors are defined, the current \code{\link{palette}} will be used.
#'
#' The time axis will be automatically labeled depending upon the length of the dendrometer series. If \code{period} is specified using a \code{numeric}, DOY values are displayed on the x-axis. In case a \code{character} of two time stamps is provided, axis labeling will be as follows: if series are longer than 120 days, years and months will be shown. If the length is between 30 and 120 days, months and days, and below 30 days, months, days and hours are displayed.
#'
#' @author Marko Smiljanic
#'
#' @references Deslauriers, A., Rossi, S., Turcotte, A., Morin, H. and Krause, C. (2011) A three-step procedure in SAS to analyze the time series from automatic dendrometers. \emph{Dendrochronologia} 29: 151-161.
#'
#' @examples
#' data(dmCD)
#' dm.phase <- phase_def(dmCD)
#' phase_plot(dmCD, dm.phase)
#'
#' # zoom in on the dendrometer series
#' phase_plot(dmCD, dm.phase, period = c(133, 142))
#'
#' # customization options
#' phase_plot(dmCD, dm.phase, period = c("2008-05-12", "2008-05-22"),
#'            colPhases = c("green", "cyan", "orange"),
#'            pch = 4, main = "Dendrometer", ylab = "Values")
#'
#' \dontrun{
#'
#' # specific sensors may be selected as follows:
#' data(dmED)
#' dm.gpf <- fill_gaps(dmED)
#' dm.phase <- phase_def(dm.gpf)
#' phase_plot(dm.gpf, dm.phase, sensor = 1)
#' phase_plot(dm.gpf, dm.phase, sensor = c(2,1))
#' phase_plot(dm.gpf, dm.phase, sensor = "Beech03")
#' phase_plot(dm.gpf, dm.phase, sensor = c("Beech03", "Beech04"))
#' }
#'
#' @import graphics
#' @import grDevices
#'
#' @export phase_plot
#'
phase_plot <- function(dm.gpf, dm.phase, sensor = NULL, period = NULL, colPhases = NULL, ...)
{
  nm1 <- deparse(substitute(dm.gpf))
  nm2 <- deparse(substitute(dm.phase))

  if(!is.dendro(dm.gpf)) {
    stop("invalid dendrometer data")
  }

  resolution <- dendro.resolution(dm.gpf)
  timestamps <- as.POSIXct(row.names(dm.gpf), tz = "GMT")

  if(!is.null(sensor)) {
    if(is.numeric(sensor)) {
      if(any(sensor > ncol(dm.gpf))) {
        sensor <- sensor[which(sensor <= ncol(dm.gpf))]
        warning("some sensors are undefined and were therefore omitted")
      }
      dm.gpf <- dm.gpf[, sensor, drop = FALSE]
      dm.phase <- dm.phase[, sensor, drop = FALSE]
    }
    else if(is.character(sensor)) {
      if(any(!(sensor %in% names(dm.gpf)))) {
        sensor <- sensor[which(sensor %in% names(dm.gpf))]
        warning("some sensors are undefined and were therefore omitted")
      }
      dm.gpf <- dm.gpf[, sensor, drop = FALSE]
      dm.phase <- dm.phase[, sensor, drop = FALSE]
    }
    else {
      warning("incorrect format of 'sensor'")
    }
  }

  if(!is.null(period) && length(period) != 2) {
    warning("'period' should have two members; the default (NULL) is used")
    period <- NULL
  }
  else if(is.numeric(period)){
    timestamp.labs <- format(timestamps, format = "%j", tz = "GMT")
    period.min <- min(period)
    period.max <- max(period)
    period.subset <- which(timestamp.labs >= period.min & timestamp.labs <= period.max)
    dm.gpf <- dm.gpf[period.subset, , drop = FALSE]
    dm.phase <- dm.phase[period.subset, , drop = FALSE]
    timestamps <- as.POSIXct(row.names(dm.gpf), tz = "GMT")
    axis.labs <- "DOY"
  }
  else if(is.character(period)) {
    options(show.error.messages = FALSE)
    period <- try(as.POSIXct(period, tz="GMT"), silent = TRUE)
    if("try-error" %in% class(period)) {
      period <- NULL
      warning("incorrect format for 'period'")
    }
    else {
      period.min <- min(period)
      period.max <- max(period)
      period.subset <- which(timestamps >= period.min & timestamps <= period.max)
      dm.gpf <- dm.gpf[period.subset, , drop = FALSE]
      dm.phase <- dm.phase[period.subset, , drop = FALSE]
      timestamps <- as.POSIXct(row.names(dm.gpf), tz = "GMT")
    }
    axis.labs <- "TIME"
  }
  else if(is.null(period)) {
    axis.labs <- "TIME"
  }

  if(length(colPhases) != 3) {
    colPhases <- NULL
  }

  dots <- list(...)
  plotArgs <- dots
  plotArgs$type <- "n"
  if(is.null(plotArgs$xaxt)) {
    plotArgs$xaxt <- "n"
  }
  else if(plotArgs$xaxt == "y") {
    plotArgs$xaxt <- NULL
  }
  else {
    plotArgs$xaxt <- "n"
  }
  if(is.null(plotArgs$las)) {
    plotArgs$las <- 1
  }

  if(class(dm.gpf) == "numeric") {
    plotArgs$x <- c(1, length(dm.gpf))
    plotArgs$y <- c(min(dm.gpf, na.rm = TRUE), max(dm.gpf, na.rm = TRUE))
  }
  else if(class(dm.gpf) == "data.frame" || class(dm.gpf) == "matrix") {
    plotArgs$x <- c(1, nrow(dm.gpf))
    plotArgs$y <- c(min(dm.gpf[1:length(dm.gpf)], na.rm = TRUE), max(dm.gpf[1:length(dm.gpf)], na.rm = TRUE))
  }
  else {
    stop(paste("'", nm1, "' must be either a numeric or data.frame", sep = ""))
  }
  if(is.null(plotArgs$ylab)) {
    plotArgs$ylab <- "Stem-size variation"
  }
  plotArgs$xlab <- ""
  do.call(plot,plotArgs)

  plotPhase <- function(den = dm.gpf, phs = dm.phase, dot = dots, colPh = colPhases) {
    pointArgs <- dot
    pointArgs$x <- 1:length(phs)
    pointArgs$y <- den[1:length(phs)]
    pointArgs$type <- "p"
    pointArgs$col <- phs
    if(!is.null(colPh)) {
      pointArgs$col <- colPh[phs]
    }
    do.call(points, pointArgs)
  }

  if(class(dm.gpf) == "numeric") {
    if(length(dm.gpf) != length(dm.phase)) {
      stop(paste("'", nm1, "' and '", nm2, "' must be of the same length", sep = ""))
    }
    lines(dm.gpf, ...)
    plotPhase()
  }
  else if(class(dm.gpf) == "data.frame" || class(dm.gpf) == "matrix") {
    dm.gpf <- as.data.frame(dm.gpf)
    if(nrow(dm.gpf) != nrow(dm.phase) || ncol(dm.gpf) != ncol(dm.phase)) {
      stop(paste("the dimensions of '", nm1, "' and '", nm2, "' must be the same", sep = ""))
    }
    for(d in 1:length(dm.gpf)) {
      dend <- dm.gpf[[d]]
      phas <- dm.phase[[d]]
      lineArgs <- dots
      lineArgs$x <- dend
      if(is.null(lineArgs$col)) {
        lineArgs$col <- d
        noColDef <- TRUE
      }
      do.call(lines,lineArgs)
      plotPhase(dend, phas)
    }
  }

  t <- axis(1, labels = FALSE)

  if(length(timestamps) < (max(t)+1)) {
    tt <- 1:(max(t)+1)
    tt[1:length(timestamps)] <- timestamps
    for(i in (length(timestamps)+1):length(tt)) {
      tt[i] <- timestamps[length(timestamps)]+(i-length(timestamps))*resolution
    }
    timestamps <- as.POSIXct(tt, origin = "1970-01-01", tz = "GMT")
  }

  axisArgs <- dots
  timestamps <- timestamps[t+1]
  timestamps <- timestamps-1*resolution
  if(axis.labs == "TIME") {
    timerange <- max(timestamps)-min(timestamps)
    if(units(timerange) != "days") {
      units(timerange) <- "days"
    }
    if(timerange > 120) {
      timestamps <- substr(timestamps, start = 1, stop = 7)
      timelabel = "[Year-Month]"
    }
    else if(timerange > 30) {
      timestamps <- substr(timestamps, start = 6, stop = 10)
      timelabel = "[Month-Day]"
    }
    else {
      timestamps <- substr(timestamps, start = 6, stop = 13)
      timelabel = "[Month-Day Hour]"
    }
    titleArgs <- dots
    if(is.null(titleArgs$xlab)) {
      titleArgs$xlab <- as.expression(bquote('Time'[.(timelabel)]))
    }
    do.call(title, titleArgs)
    if(!is.null(axisArgs$cex.axis)) {
      axisArgs$cex.axis = axisArgs$cex.axis * 0.75
    }
    else {
      axisArgs$cex.axis = par("cex.axis") * 0.75
    }

  }
  else if(axis.labs == "DOY") {
    timestamps <- format(timestamps, format = "%j", tz = "GMT")
    titleArgs <- dots
    if(is.null(titleArgs$xlab)) {
      titleArgs$xlab <- as.expression(bquote('Day of year'))
    }
    do.call(title, titleArgs)

  }

  axisArgs$side <- 1
  axisArgs$at <- t
  axisArgs$labels <- timestamps
  axisArgs$col <- NULL
  do.call(axis, axisArgs)
}
