#' Solar angles
#'
#' This function returns the solar angles for a given time and location.
#'
#' @param time POSIXct Time, any valid time zone (TZ) is allowed, default is
#'   current time
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees).
#' @param lon numeric Vector of longitudes (degrees)
#' @param lat numeric Vector of latitudes (degrees)
#' @param use_refraction logical Flag indicating whether to correct for
#'   fraction in the atmosphere
#'
#' @return A list with components time in same TZ as input, azimuth, elevation,
#'   diameter, and distance.
#'
#' @family astronomy related functions
#'
#'
#' @export
#' @examples
#' require(lubridate)
#' sun_angles()
#' sun_angles(ymd_hms("2014-09-23 12:00:00"))
#' sun_angles(ymd_hms("2014-09-23 12:00:00"), lat=60, lon=0)
#'
sun_angles <- function(time = lubridate::now(),
                       geocode = NULL,
                       lon = 0, lat = 0,
                       use_refraction = FALSE)
{
  stopifnot(lubridate::is.POSIXct(time))
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  tz <- lubridate::tz(time)
  t <- lubridate::with_tz(time, "UTC")
  nt <- length(t)

  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  nlon <- length(lon)
  nlat <- length(lat)
  if (nlon != nlat)
    stop("lengths of longitude and latitude must match")
  if (nlon == 1) {
    lon <- rep(lon, nt)
    lat <- rep(lat, nt)
  }
  else {
    if (nt != nlon)
      stop("lengths of t, latitude and longitude must match, unless last two are of length 1")
  }
  year <- lubridate::year(t)
  if (any(year < 1950) || any(year > 2050))
    stop("year=", year, " is outside acceptable range")
  day <- lubridate::yday(t)
  if (any(day < 1) || any(day > 366))
    stop("day is not in range 1 to 366")
  hour <- lubridate::hour(t) + lubridate::minute(t) / 60 + lubridate::second(t) / 3600
  if (any(hour < -13) || any(hour > 36))
    stop("hour outside range -13 to 36")
  if (any(lat < -90)) {
    warning("latitude(s) trimmed to range -90 to 90")
    lat[lat < -90] <- -90
  }
  if (any(lat > 90)) {
    warning("latitude(s) trimmed to range -90 to 90")
    lat[lat > 90] <- 90
  }
  if (any(lon < -180)) {
    warning("longitude(s) trimmed to range -180 to 180")
    lon[lon < -180] <- -180
  }
  if (any(lon > 180)) {
    warning("longitude(s) trimmed to range -180 to 180")
    lon[lon > 180] <- 180
  }
  delta <- year - 1949
  leap <- delta %/% 4
  jd <- 32916.5 + (delta * 365 + leap + day) + hour / 24
  jd <- jd + ifelse(0 == (year %% 100) & 0 != (year %% 400), 1,
                    0)
  time <- jd - 51545
  mnlong <- 280.46 + 0.9856474 * time
  mnlong <- mnlong %% 360
  mnlong <- mnlong + ifelse(mnlong < 0, 360, 0)
  mnanom <- 357.528 + 0.9856003 * time
  mnanom <- mnanom %% 360
  mnanom <- mnanom + ifelse(mnanom < 0, 360, 0)
  rpd <- pi/180
  mnanom <- mnanom * rpd
  eclong <- mnlong + 1.915 * sin(mnanom) + 0.02 * sin(2 * mnanom)
  eclong <- eclong %% 360
  eclong <- eclong + ifelse(eclong < 0, 360, 0)
  oblqec <- 23.439 - 4e-07 * time
  eclong <- eclong * rpd
  oblqec <- oblqec * rpd
  num <- cos(oblqec) * sin(eclong)
  den <- cos(eclong)
  ra <- atan(num/den)
  ra <- ra + ifelse(den < 0, pi, ifelse(num < 0, 2 * pi, 0))
  dec <- asin(sin(oblqec) * sin(eclong))
  gmst <- 6.697375 + 0.0657098242 * time + hour
  gmst <- gmst %% 24
  gmst <- gmst + ifelse(gmst < 0, 24, 0)
  lmst <- gmst + lon/15
  lmst <- lmst %% 24
  lmst <- lmst + ifelse(lmst < 0, 24, 0)
  lmst <- lmst * 15 * rpd
  ha <- lmst - ra
  ha <- ha + ifelse(ha < (-pi), 2 * pi, 0)
  ha <- ha - ifelse(ha > pi, 2 * pi, 0)
  el <- asin(sin(dec) * sin(lat * rpd) + cos(dec) * cos(lat *
                                                          rpd) * cos(ha))
  az <- asin(-cos(dec) * sin(ha)/cos(el))
  az <- ifelse(sin(dec) - sin(el) * sin(lat * rpd) > 0,
               ifelse(sin(az) < 0, az + 2 * pi, az), pi - az)
  el <- el/rpd
  az <- az/rpd
  if (use_refraction) {
    refrac <- ifelse(el >= 19.225, 0.00452 * 3.51823/tan(el * rpd),
                     ifelse(el > (-0.766) & el < 19.225,
                            3.51823 * (0.1594 + el * (0.0196 + 2e-05 * el))/(1 + el * (0.505 + 0.0845 * el)),
                            0))
    el <- el + refrac
  }
  soldst <- 1.00014 - 0.01671 * cos(mnanom) - 0.00014 * cos(2 * mnanom)
  soldia <- 0.5332/soldst
  if (any(el < (-90)) || any(el > 90))
    stop("output el out of range")
  if (any(az < 0) || any(az > 360))
    stop("output az out of range")
  return(list(time = lubridate::with_tz(t, tz),
              azimuth = az,
              elevation = el,
              diameter = soldia,
              distance = soldst))
}


#' Times for sun positions
#'
#' Functions for calculating the timing of solar positions by means of function
#' \code{sun_angles}, given geographical coordinates and dates. They can be also
#' used to find the time for an arbitrary solar elevation between 90 and -90
#' degrees by supplying "twilight" angle(s) as argument.
#'
#' @param date array of POSIXct times or Date objects, any valid TZ is allowed,
#'   default is current date
#' @param tz character string incading time zone to be used in output, default
#'   is system time zone
#' @param geocode data frame with variables lon and lat as numeric values
#'   (degrees).
#' @param lon numeric array of longitudes (degrees)
#' @param lat numeric array of latitudes (degrees)
#' @param twilight character string, one of "none", "civil", "nautical",
#'   "astronomical", or a \code{numeric} vector of length one, or two, giving
#'   solar elevation angle(s) in degrees (negative if below the horizon).
#' @param unit.out charater string, One of "date", "hour", "minute", or "second".
#'
#' @return \code{day_night} returns a list with fields sunrise time, sunset
#'   time, day length, night length. Each element of the list is a vector of the
#'   same length as the argument supplied for date.
#'
#' @note If twilight is a numeric vector of length two, the element with index 1
#'   is used for sunrise and that with index 2 for sunset.
#'
#' @family astronomy related functions
#'
#' @name day_night
#' @export
#' @examples
#' library(lubridate)
#' day_length()
#' day_length(ymd("2015-05-30"), lat = 60, lon = 25)
#' day_length(ymd("2014-12-30"), lat = 60, lon = 25)
#' day_length(ymd("2015-05-30"), lat = 60, lon = 25, twilight = "civil")
#' sunrise_time(ymd("2015-05-30"), lat = 60, lon = 25, tz = "EET")
#' day_night(ymd("2015-05-30"), lat = 60, lon = 25, twilight = "civil")
#'
day_night <- function(date = lubridate::today(),
                      tz = "UTC",
                      geocode = NULL, lon = 0, lat = 0,
                      twilight = "none",
                      unit.out = "date") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }

  list(day         = as.Date(date),
       sunrise     = sunrise_time(date = date, tz = tz,
                                  lon = lon, lat = lat,
                                  twilight = twilight,
                                  unit.out = unit.out),
       noon        = noon_time(date = date, tz = tz,
                               lon = lon, lat = lat,
                               unit.out = unit.out),
       sunset      = sunset_time(date = date, tz = tz,
                                 lon = lon, lat = lat,
                                 twilight = twilight,
                                 unit.out = unit.out),
       daylength   = day_length(date = date, tz = tz,
                                lon = lon, lat = lat,
                                twilight = twilight,
                                unit.out = unit.out),
       nightlength = night_length(date = date, tz = tz,
                                  lon = lon, lat = lat,
                                  twilight = twilight,
                                  unit.out = unit.out) )
}

#' twilight argument check and conversion
#'
#' @return numeric  Solar elevation angle at sunrise or sunset
#' @keywords internal
twilight2angle <- function(twilight) {
  if (!is.numeric(twilight)) {
    if (twilight == "none") {
      twilight_angle <- c(0, 0)
    } else if (twilight == "civil") {
      twilight_angle <- c(-6, -6)
    } else if (twilight == "nautical") {
      twilight_angle <- c(-12, -12)
    } else if (twilight == "astronomical") {
      twilight_angle <- c(-18, -18)
    } else {
      twilight_angle <- c(NA, NA)
    }
  } else {
    if (length(twilight) == 1) {
      twilight_angle <- rep(twilight, 2)
    } else if (length(twilight) == 2) {
      twilight_angle <- twilight
    } else {
      twilight_angle <- c(NA, NA)
    }
    twilight_angle <- ifelse(twilight_angle < 90, twilight_angle, NA)
    twilight_angle <- ifelse(twilight_angle > -90, twilight_angle, NA)
  }
  if (any(is.na(twilight_angle))) {
    stop("Unrecognized argument value for 'twilight': ", twilight)
  }
  twilight_angle
}

#' date argument check and conversion
#'
#' @return numeric representtaion of the date
#' @keywords internal
date2seconds <- function(t, tz) {
  if (!lubridate::is.POSIXct(t)) {
    if (lubridate::is.instant(t)) {
      t <- as.POSIXct(t, tz = "UTC")
    } else {
      warning("t is not a valid time or date")
    }
  }
  t <- as.POSIXct(t, tz = tz)
  lubridate::hour(t) <- 0
  lubridate::minute(t) <- 0
  lubridate::second(t) <- 0
  as.numeric(t, tz = tz)
}

#' time argument check and conversion
#'
#' @return numeric representtaion of the date
#' @keywords internal
time2seconds <- function(t, tz) {
  if (!lubridate::is.POSIXct(t)) {
    if (lubridate::is.instant(t)) {
      t <- as.POSIXct(t, tz = "UTC")
    } else {
      warning("t is not a valid time or date")
    }
  }
  t <- as.POSIXct(t, tz = tz)
  as.numeric(t, tz = tz)
}

#' function to be numerically minimized
#'
#' @return an elevation angle delta
#' @keywords internal
altitude <- function(x, lon, lat, twlght_angl = 0){
  t_temp <- as.POSIXct(x, origin = lubridate::origin, tz = "UTC")
  return(sun_angles(t_temp,
                    lon = lon,
                    lat = lat)$elevation - twlght_angl)
}

#' @rdname day_night
#' @export
#' @return \code{noon_time}, \code{sunrise_time} and \code{sunset_time} return a
#'   vector of POSIXct times
noon_time <- function(date = lubridate::today(), tz = "UTC",
                      geocode = NULL, lon = 0, lat = 0,
                      twilight = NA, unit.out = "date") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  date_num <- sapply(date, date2seconds, tz = tz)
  times <- numeric()
  twlght_angl <- 0
  for (t_num in date_num) {
    noon <- try(
      stats::optimize(f = altitude, interval = c(t_num + 7200, t_num + 86400 - 7200),
               lon = lon, lat = lat,
               twlght_angl = twlght_angl,
               maximum = TRUE)$maximum
    )
    if (inherits(noon, "try-error")) {
      noon <- NA
    }
    times <- c(times, noon)
  }
  times <- as.POSIXct(times, tz = tz, origin = lubridate::origin)
  if (unit.out != "date") {
    times <- sapply(times, date2tod, unit.out = unit.out)
  }
  times
}

#' @rdname day_night
#'
#' @export
sunrise_time <- function(date = lubridate::today(), tz = "UTC",
                         geocode = NULL, lon = 0, lat = 0,
                         twilight = "none", unit.out = "date") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  noon <- noon_time(date = date, tz = tz,
                    lon = lon, lat = lat)
  noon_num <- sapply(noon, time2seconds, tz = tz)
  times <- numeric()
  twlght_angl <- twilight2angle(twilight)[1]
  for (t_num in noon_num) {
     rise <- try(
      stats::uniroot(f = altitude,
              lon = lon, lat = lat, twlght_angl = twlght_angl,
              lower = t_num - 86400/2, upper = t_num)$root,
      silent = TRUE)
    if (inherits(rise, "try-error")) {
      rise <- NA # never
    }
    times <- c(times, rise)
  }
  times <- as.POSIXct(times, tz = tz, origin = lubridate::origin)
  if (unit.out != "date") {
    times <- sapply(times, date2tod, unit.out = unit.out)
  }
  times
}

#' @rdname day_night
#' @export
#'
#' @note \code{night_length} returns the length of night-time conditions in one
#'   day (00:00:00 to 23:59:59), rather than the length of the night between two
#'   consequtive days.
sunset_time <- function(date = lubridate::today(), tz = "UTC",
                        geocode = NULL, lon = 0, lat = 0,
                        twilight = "none", unit.out = "date") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  noon <- noon_time(date = date, tz = tz,
                    lon = lon, lat = lat)
  noon_num <- sapply(noon, time2seconds, tz = tz)
  times <- numeric()
  twlght_angl <- twilight2angle(twilight)[2]
  for (t_num in noon_num) {
    set <- try(
      stats::uniroot(f = altitude,
              lon = lon, lat = lat,
              twlght_angl = twlght_angl,
              lower = t_num, upper = t_num + 86400 / 2)$root,
      silent = TRUE)
    if (inherits(set, "try-error")) {
      set <- NA # never
    }
    times <- c(times, set)
  }
  times <- as.POSIXct(times, tz = tz, origin = lubridate::origin)
  if (unit.out != "date") {
    times <- sapply(times, date2tod, unit.out = unit.out)
  }
  times
}

#' @rdname day_night
#'
#' @export
#' @return \code{day_length} and \code{night_length} return numeric a vector
#'   giving the length in hours
day_length <- function(date = lubridate::today(), tz = "UTC",
                       geocode = NULL, lon = 0, lat = 0,
                       twilight = "none", unit.out = "hour") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  noon <- noon_time(date = date, tz = tz,
                    lon = lon, lat = lat)
  rise_time <- sunrise_time(date = date, tz = tz,
                            lon = lon, lat = lat,
                            twilight = twilight, unit.out = "date")
  set_time <- sunset_time(date = date, tz = tz,
                          lon = lon, lat = lat,
                          twilight = twilight, unit.out = "date")
  hours <- ifelse(is.na(rise_time) | is.na(set_time),
         ifelse(altitude(noon, lon = lon, lat = lat) > 0, 24, 0),
         set_time - rise_time)
  switch(unit.out,
         "date" = hours,
         "hour" = hours,
         "minute" = hours * 60,
         "second" = hours * 3600)
}

#' @rdname day_night
#'
#' @export
night_length <- function(date = lubridate::today(), tz = "UTC",
                         geocode = NULL, lon = 0, lat = 0,
                         twilight = "none", unit.out = "hour") {
  stopifnot(is.null(geocode) || is.data.frame(geocode))
  if (!is.null(geocode)) {
    lon <- geocode[["lon"]]
    lat <- geocode[["lat"]]
    geocode <- NULL
  }
  hours <- 24 - day_length(date = date, tz = tz,
                           lon = lon, lat = lat,
                           twilight = twilight, unit.out = "hour")
  switch(unit.out,
         "date" = hours,
         "hour" = hours,
         "minute" = hours * 60,
         "second" = hours * 3600)
}

#' Convert date to time-of-day in hours
#'
#' @param date a date object accepted by lubridate functions
#'
#' @keywords internal
date2tod <- function(date, unit.out) {
  if (unit.out == "hour") {
    lubridate::hour(date) + lubridate::minute(date) / 60 + lubridate::second(date) / 3600
  } else if (unit.out == "minute") {
    lubridate::hour(date) * 60 + lubridate::minute(date) + lubridate::second(date) / 60
  } else if (unit.out == "second") {
    lubridate::hour(date) * 3600 + lubridate::minute(date) * 60 + lubridate::second(date)
  } else {
    stop("Unrecognized 'unit.out': ", unit.out)
  }
}

