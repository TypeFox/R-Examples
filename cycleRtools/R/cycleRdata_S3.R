#' Summary method for cycleRdata class.
#'
#' Relevant summary metrics for cycling data (method for class
#' \code{"cycleRdata"}).
#'
#' @param object object for which a summary is desired.
#' @param sRPE optional; session Rating of Percieved Exertion (value between 1
#'   and 10; Foster 1998).
#' @param CP optional; Critical Power value (Watts).
#' @param .smoothpwr character string; column name of smoothed power values.
#'   Used for xP metric.
#' @param ... further arguments passed to or from other methods.
#'
#' @return a list object of class \code{"cyclesummary"}, which has an associated
#'   print method.
#'
#' @examples
#' data(intervaldata)
#' summary(intervaldata)
#'
#' @references Foster C. Monitoring training in athletes with reference to
#' overtraining syndrome. Medicine & Science in Sports & Exercise 30: 1164-1168,
#' 1998.
#'
#' @export
summary.cycleRdata <- function(object, sRPE = attr(object, "sRPE"), CP = attr(object, "CP"),
                               .smoothpwr = "power.smooth.W", ...) {
  if (!is.null(sRPE) && (sRPE > 10 || sRPE < 0)) sRPE <- NULL
  # Functions.
  defaultNA <- function(need, expr) if (need) eval(expr) else NA
  xP        <- function(x) pwr_tf(x, 4)
  # Summary functions with na.rm = TRUE default.
  my        <- list(mean = function(x) mean(x, na.rm = TRUE),
                    sd   = function(x) sd(x, na.rm = TRUE),
                    max  = function(x) max(x, na.rm = TRUE),
                    sum  = function(x) sum(x, na.rm = TRUE),
                    tapl = function(x, f) tapply(object[[x]], object[["lap"]], f),
                    tapi = function(x, f) tapply(object[[x]], object[["interval"]], f))
  climbin   <- object$delta.elev > 0  # For brevity.
  tests     <- list(dist = "distance.km", spd = "speed.kmh", ele = "elevation.m",
                    work = "work.kJ", pwr = "power.W", smth = .smoothpwr)
  tests     <- lapply(tests, function(x) all(!is.na(object[, x])))
  metrics   <- list()
  # Date-time metrics ----------------------------------------------------------
  metrics$date_time <- data.frame(
    ride.date        = format(object$timestamp[[1]], format = "%d %b %Y"),
    ride.start.time  = format(object$timestamp[[1]], format = "%H:%M:%S")
  )
  # Overall metrics ------------------------------------------------------------
  needs   <- with(tests, list(  # Assign names here.
    ride.time.min    = TRUE,
    elapsed.time.min = TRUE,
    distance.km      = dist,
    work.kJ          = work,
    norm.work.kJ     = pwr,
    climbing.m       = ele,
    sRPE.score       = !is.null(sRPE),
    TRIMP.score      = !is.null(CP) && pwr,
    speed.mean.kmh   = spd,
    speed.sd.kmh     = spd,
    power.mean.W     = pwr,
    power.sd.W       = pwr,
    xPower.W         = smth,
    power.max.W      = pwr
  ))
  exprs   <- alist(
    ride_time(object$timer.s) / 60,
    my$max(object$timer.s) / 60,
    my$max(object$distance.km),
    my$max(object$work.kJ),
    (ride_time(object$timer.s) * xP(object[[.smoothpwr]])) / 1000,
    my$sum(object[climbin, "delta.elev"]),
    sRPE * (ride_time(object$timer.s) / 60),
    pwr_TRIMP(object, CP),
    my$mean(object$speed.kmh),
    my$sd(object$speed.kmh),
    my$mean(object$power.W),
    my$sd(object$power.W),
    xP(object[[.smoothpwr]]),
    my$max(object$power.W)
  )
  metrics$overall <-
    data.frame(mapply(defaultNA, needs, exprs, SIMPLIFY = FALSE))
  # Power-time profile ---------------------------------------------------------
  if (tests$pwr) {
    wins       <- c(1, 2, 5, 10, 15, 20) * 60
    metrics$Pt <- mmv(object, "power.W", wins)[1, ]
  }
  # Lap metrics ----------------------------------------------------------------
  if (length(unique(object$lap)) > 1) {
    needs <- with(tests, list(
      lap.no             = TRUE,
      lap.time.min       = TRUE,
      lap.distance.km    = dist,
      lap.work.kJ        = work,
      lap.climbing.m     = ele,
      lap.speed.mean.kmh = spd,
      lap.speed.sd.kmh   = spd,
      lap.power.mean.W   = pwr,
      lap.power.sd.W     = pwr,
      lap.xPower.W       = smth,
      lap.power.max.W    = pwr
    ))
    exprs <- alist(
      as.numeric(unique(object$lap)),
      my$tapl("timer.s", ride_time) / 60,
      my$tapl("distance.km", my$max),
      my$tapl("work.kJ", my$max),
      tapply(object[climbin, "delta.elev"], object[climbin, "lap"], my$sum),
      my$tapl("speed.kmh", my$mean),
      my$tapl("speed.kmh", my$sd),
      my$tapl("power.W", my$mean),
      my$tapl("power.W", my$sd),
      my$tapl(.smoothpwr, xP),
      my$tapl("power.W", my$max)
    )
    metrics$laps <-
      data.frame(mapply(defaultNA, needs, exprs, SIMPLIFY = FALSE))
    metrics$laps$lap.distance.km <-
      with(metrics$laps, c(lap.distance.km[1], Diff(lap.distance.km)))
    metrics$laps$lap.work.kJ <-
      with(metrics$laps, c(lap.work.kJ[1], Diff(lap.work.kJ)))
  }
  # Interval metrics -----------------------------------------------------------
  if (!is.null(object$interval)) {
    needs <- with(tests, list(
      interval.no             = TRUE,
      interval.time.min       = TRUE,
      interval.distance.km    = dist,
      interval.work.kJ        = work,
      interval.climbing.m     = ele,
      interval.speed.mean.kmh = spd,
      interval.speed.sd.kmh   = spd,
      interval.power.mean.W   = pwr,
      interval.power.sd.W     = pwr,
      interval.xPower.W       = smth,
      interval.power.max.W    = pwr
    ))
    exprs <- alist(
      as.numeric(unique(object$interval)),
      my$tapi("timer.s", ride_time) / 60,
      my$tapi("distance.km", my$max),
      my$tapi("work.kJ", my$max),
      tapply(object[climbin, "delta.elev"], object[climbin, "interval"], my$sum),
      my$tapi("speed.kmh", my$mean),
      my$tapi("speed.kmh", my$sd),
      my$tapi("power.W", my$mean),
      my$tapi("power.W", my$sd),
      my$tapi(.smoothpwr, xP),
      my$tapi("power.W", my$max)
    )
    metrics$intervals <-
      data.frame(mapply(defaultNA, needs, exprs, SIMPLIFY = FALSE))
    metrics$intervals$interval.distance.km <-
      with(metrics$intervals, c(interval.distance.km[1], Diff(interval.distance.km)))
    metrics$intervals$interval.work.kJ <-
      with(metrics$intervals, c(interval.work.kJ[1], Diff(interval.work.kJ)))
  }
  class(metrics) <- c("cyclesummary", class(metrics))
  metrics
}
#' @export
print.cyclesummary <- function(x, ...) {
  cat("\n")
  cat(paste("Date:      ", x$date_time$ride.date, "\n"))
  cat(paste("Ride start:", x$date_time$ride.start.time, "\n"))

  cat("\n**SUMMARY**\n")
  table <- t(x$overall); colnames(table) <- " "
  print(table, digits = 2)

  if (!is.null(x$Pt)) {
    cat("\n**POWER-TIME PROFILE**\n\n")
    table <- matrix(x$Pt, ncol = 1, dimnames = list(
      paste(as.numeric(names(x$Pt)) / 60, "min"), "Watts"))
    print(table, digits = 2)
  }

  if (!is.null(x$laps)) {
    cat("\n**LAPS**\n\n")
    table <- x$laps; colnames(table) <- c("#", gsub("lap.", "", colnames(table))[-1])
    print(table, digits = 2, row.names = FALSE)
  }

  if (!is.null(x$intervals)) {
    cat("\n**INTERVALS**\n\n")
    table <- x$intervals; colnames(table) <- c("#", gsub("interval.", "", colnames(table))[-1])
    print(table, digits = 2, row.names = FALSE)
  }
}
#' Plot cycling data.
#'
#' Generate plots to effectively summarise a cycling dataset.
#'
#' The \code{y} argument describes plot options such that: \enumerate{\item
#' plots W' balance (kJ). \item plots power data (W). \item plots an elevation
#' profile (m).} These options can be combined to produce a stack of plots as
#' desired.
#'
#' @param x a \code{"cycleRdata"} object produced by \code{read*()}.
#' @param y numeric; plots to be created (see details).
#' @param xvar character; name of the column to be plotted as the xvariable.
#' @param xlab character; x axis label for bottom plot.
#' @param xlim given in terms of \code{x}.
#' @param CP a value for critical power annotation.
#' @param laps logical; should laps be seperately coloured?
#' @param breaks logical; should plot lines be broken when stationary? Will only
#'   show when \code{xvar} represents time values.
#' @param ... graphical parameters, and/or arguments to be passed to or from
#'   other methods.
#'
#' @return a variable number of plots.
#'
#' @examples \dontrun{
#' data(ridedata)
#'
#' plot(ridedata, xvar = "timer.min")
#' plot(ridedata, xvar = "distance.km")
#'
#' ## With only two plots.
#' plot(ridedata, y = c(2, 1))
#'
#' ## Using xlim, note that title metrics adjust.
#' plot(ridedata, xvar = "timer.min", xlim = c(100, 150))
#'
#' ## Lap colouring.
#' data(intervaldata)
#' plot(intervaldata, xvar = "timer.min", laps = TRUE)
#' }
#'
#' @export
plot.cycleRdata <- function(x, y = 1:3, xvar = "timer.s", xlab = NULL, xlim = NULL,
                            CP = attr(x, "CP"), laps = FALSE, breaks = TRUE, ...) {
  ## See R/special_plots.R ##
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar[-grep("pin", names(opar))]))

  if (breaks)
    x <- expand_stops(x, keep_attr = TRUE)

  if (!is.null(xlim))  # Subset x according to xlim arg (for titles).
    x <- x[x[, xvar] %btwn% xlim, ]

  if (length(y) == 1) {  # Single plot.
    par(...)
    switch(y,
           Wbal_plot(x, xvar, CP),
           pwr_plot(x, xvar, laps, CP),
           elev_plot(x, xvar))

  } else {  # Stack of plots.

    dev <- try(getOption("device"))
    if (!is.function(dev) && dev == "RStudioGD") {
      if (!identical(y, (y <- test_y(y)))) warning(paste(
        "Graphics device not big enough for requested plots."), call. = FALSE)
    }

    par(mfrow = c(length(y), 1), ...)
    for (i in y)
      switch(i,
             Wbal_plot(x, xvar, CP),
             pwr_plot(x, xvar, laps, CP),
             elev_plot(x, xvar))
  }

  # Generate bottom most xlab.
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = 2.5)
  else
    switch(xvar,
           "timer.s"     = mtext("Time (sec)",    side = 1, line = 2.5),
           "timer.min"   = mtext("Time (min)",    side = 1, line = 2.5),
           "distance.km" = mtext("Distance (km)", side = 1, line = 2.5)
    )
}
test_y <- function(y) {
  repeat {
    par(mfrow = c(length(y), 1))
    plexpr <- substitute( # Invisible plot.
      plot(1:10, ylab = "W' expended (kJ)", xlab = "",
           main = "This is a test plot.", bty = "n", axes = FALSE,
           col = "white", col.axis = "white",
           col.lab = "white", col.main = "white"))
    test <- tryCatch(eval(plexpr), error = function(e) return("error"))
    if (is.null(test)) break else y <- y[-length(y)]
  }
  y
}
# ----------------------------------------------- #
## Test class.
#' @rdname cycleRdata
#' @export
is.cycleRdata <- function(x) {
  all(cycleRdata_fields() %in% names(x)) && any(class(x) == "cycleRdata")
}

## Coercion functions.
#' @rdname cycleRdata
#' @export
as.cycleRdata <- function(x) UseMethod("as.cycleRdata")
#' @export
as.cycleRdata.default <- function(x) {
  warning("Only data.frames can be coerced to cycleRdata class.")
  return(x)
}
#' @export
as.cycleRdata.data.frame <- function(x) {
  if (all(cycleRdata_fields() %notin% names(x)))
    warning("Required column names not detected, can't cooerce.", call. = FALSE)
  else
    class(x) <- c("cycleRdata", "data.frame")
  x
}

## Protect column names.
#' @export
`names<-.cycleRdata` <- function(x, value) {
  class(x) <- "data.frame"
  names(x) <- value
  if (all(cycleRdata_fields() %notin% names(x)))
    warning("Improper names found, removing 'cycleRdata' class attribute.", call. = FALSE)
  else
    class(x) <- c("cycleRdata", "data.frame")
  x
}
