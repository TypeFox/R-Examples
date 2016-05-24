# Column names of cycleRdata objects, for the reference.
cycleRdata_fields <- function() {
  c("timer.s", "timer.min", "timestamp", "delta.t",    # Time.
    "lat", "lon", "distance.km", "speed.kmh",          # GPS.
    "elevation.m", "delta.elev", "VAM",                # Altitude.
    "power.W", "power.smooth.W", "work.kJ", "Wexp.kJ", # Power.
    "cadence.rpm", "hr.bpm", "lap")                    # Other.
}
eval_check <- function(expr, envir) {  # Handles missing columns.
  out <- try(eval(parse(text = expr), envir), silent = TRUE)
  if (inherits(out, "try-error")) out <- rep_len(NA, nrow(envir))
  out
}
timediff <- function(x) {
  c(0, Diff(x - x[[1]]))
}
## Captured expressions are quoted and parsed so that they don't raise
## errors on R CMD check.
# ------------------------------------------------------------------------------
format_fit <- function(data) {
  colnames(data) <- gsub("^record.", "", colnames(data))  # For brevity.
  exprs <- alist(
    timer.s        = "timestamp.s. - timestamp.s.[[1]]",
    timer.min      = "(timestamp.s. - timestamp.s.[[1]]) / 60",
    timestamp      = "as.POSIXct(timestamp.s., tz = 'UTC', origin = '1989-12-31')",
    delta.t        = "timediff(timestamp.s.)",
    lat            = "position_lat.semicircles. * (180 / (2 ^ 31))",
    lon            = "position_long.semicircles. * (180 / (2 ^ 31))",
    distance.km    = "distance.m. / 1000",
    speed.kmh      = "(speed.m.s. * 60^2) / 1000",
    elevation.m    = "altitude.m.",
    delta.elev     = "c(0, Diff(altitude.m.))",
    VAM            = "NA",
    power.W        = "power.watts.",
    power.smooth.W = "NA",
    work.kJ        = "cumsum(power.watts. * timediff(timestamp.s.)) / 1000",
    Wexp.kJ        = "NA",
    cadence.rpm    = "cadence.rpm.",
    hr.bpm         = "heart_rate.bpm.",
    lap            = "lap"
  )
  out <- data.frame(lapply(exprs, eval_check, envir = data))
  # Fill in columns.
  out$VAM <- c(0, out$delta.elev[-1] / out$timer.s[-1]); out$VAM[out$VAM < 0] <- 0
  out$power.smooth.W <- rollmean_smth(as.cycleRdata(out), "power.W", 25, ema = TRUE)
  out
}
format_pwx <- function(data) {
  exprs <- alist(
    timer.s        = "timeoffset - timeoffset[[1]]",
    timer.min      = "(timeoffset - timeoffset[[1]]) / 60",
    timestamp      = "timestamp",
    delta.t        = "timediff(timeoffset)",
    lat            = "NA",
    lon            = "NA",
    distance.km    = "dist / 1000",
    speed.kmh      = "(spd * 60^2) / 1000",
    elevation.m    = "predict(smooth.spline(alt))$y",
    delta.elev     = "c(0, Diff(predict(smooth.spline(alt))$y))",
    VAM            = "NA",
    power.W        = "pwr",
    power.smooth.W = "NA",
    work.kJ        = "cumsum(pwr * timediff(timeoffset)) / 1000",
    Wexp.kJ        = "NA",
    cadence.rpm    = "cad",
    hr.bpm         = "hr",
    lap            = "1"
  )
  out <- data.frame(lapply(exprs, eval_check, envir = data))
  # Fill in columns.
  out$VAM <- c(0, out$delta.elev[-1] / out$timer.s[-1]); out$VAM[out$VAM < 0] <- 0
  out$power.smooth.W <- rollmean_smth(as.cycleRdata(out), "power.W", 25, ema = TRUE)
  out
}
format_tcx <- function(data) {
  data$Time <- as.numeric(data$Time)
  exprs <- alist(
    timer.s        = "Time - Time[[1]]",
    timer.min      = "(Time - Time[[1]]) / 60",
    timestamp      = "Time",
    delta.t        = "timediff(Time)",
    lat            = "LatitudeDegrees",
    lon            = "LongitudeDegrees",
    distance.km    = "DistanceMeters / 1000",
    speed.kmh      = "(Speed * 60^2) / 1000",
    elevation.m    = "AltitudeMeters",
    delta.elev     = "c(0, Diff(AltitudeMeters))",
    VAM            = "NA",
    power.W        = "Watts",
    power.smooth.W = "NA",
    work.kJ        = "cumsum(Watts * timediff(Time)) / 1000",
    Wexp.kJ        = "NA",
    cadence.rpm    = "Cadence",
    hr.bpm         = "HeartRateBpm",
    lap            = "1"
  )
  out <- data.frame(lapply(exprs, eval_check, envir = data))
  # Fill in columns.
  out$VAM <- c(0, out$delta.elev[-1] / out$timer.s[-1]); out$VAM[out$VAM < 0] <- 0
  out$power.smooth.W <- rollmean_smth(as.cycleRdata(out), "power.W", 25, ema = TRUE)
  out
}
format_srm <- function(data) {
  exprs <- alist(
    timer.s        = "time - time[[1]]",
    timer.min      = "(time - time[[1]]) / 60",
    timestamp      = "as.POSIXct(time, tz = '', origin = '1970-01-01')",
    delta.t        = "timediff(time)",
    lat            = "NA",
    lon            = "NA",
    distance.km    = "cumsum(((speed * 60^2) / 1000) * (timediff(time) / (60 ^ 2)))",
    speed.kmh      = "(speed * 60^2) / 1000",
    elevation.m    = "ele",
    delta.elev     = "c(0, Diff(ele))",
    VAM            = "NA",
    power.W        = "pwr",
    power.smooth.W = "NA",
    work.kJ        = "cumsum(pwr * timediff(time)) / 1000",
    Wexp.kJ        = "NA",
    cadence.rpm    = "cad",
    hr.bpm         = "hr",
    lap            = "1"
  )
  out <- data.frame(lapply(exprs, eval_check, envir = data))
  # Fill in columns.
  out$VAM <- c(0, out$delta.elev[-1] / out$timer.s[-1]); out$VAM[out$VAM < 0] <- 0
  out$power.smooth.W <- rollmean_smth(as.cycleRdata(out), "power.W", 25, ema = TRUE)
  out
}
format_GC <- function(data, filename) {
  exprs <- alist(
    timer.s        = "secs",
    timer.min      = "secs / 60",
    timestamp      = "NA",
    delta.t        = "c(0, Diff(secs))",
    lat            = "lat",
    lon            = "lon",
    distance.km    = "km",
    speed.kmh      = "kph",
    elevation.m    = "alt",
    delta.elev     = "c(0, Diff(alt))",
    VAM            = "NA",
    power.W        = "watts",
    power.smooth.W = "NA",
    work.kJ        = "cumsum(watts * c(0, Diff(secs))) / 1000",
    Wexp.kJ        = "NA",
    cadence.rpm    = "cad",
    hr.bpm         = "hr",
    lap            = "1"
  )
  out <- data.frame(lapply(exprs, eval_check, envir = data))
  # Fill in columns.
  out$VAM <- c(0, out$delta.elev[-1] / out$timer.s[-1]); out$VAM[out$VAM < 0] <- 0
  out$power.smooth.W <- rollmean_smth(as.cycleRdata(out), "power.W", 25, ema = TRUE)
  out$timestamp <- as.POSIXct(tools::file_path_sans_ext(basename(filename)),
                              format = "%Y_%m_%d_%H_%M_%S", tz = "")
  out$timestamp <- out$timestamp + out$timer.s
  out
}
