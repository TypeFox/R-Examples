#' cycleRdata class
#'
#' A class for imported ride files intended to ease integration with package
#' functionality. Produced by invoking \code{\link{read_ride}} (or equivalent)
#' with the argument \code{format = TRUE}. Fundamentally, \code{cycleRdata}
#' objects are a special type of \code{data.frame}; special in the sense that
#' column names are predefined and assumed to be present in the class'
#' associated methods. Modification of these column names will lead to errors.
#' See below for a description of the format.
#'
#' @param x an object to be tested/coerced.
#'
#' @format The columns of cycleRdata objects are structured as such:
#' \describe{
#' \item{timer.s}{an ongoing timer (seconds). Stoppages are not recorded per se,
#' but rather represented as breaks in the continuity of the timer.}
#' \item{timer.min}{as above, but in units of minutes.}
#' \item{timestamp}{"POSIXct" values, describing the actual time of day.}
#' \item{delta.t}{delta time values (seconds).}
#' \item{lat}{latitude values (degrees).}
#' \item{lng}{longitude values (degrees).}
#' \item{distance.km}{cumulative distance (kilometres).}
#' \item{speed.kmh}{speed in kilometres per hour.}
#' \item{elevation.m}{altitude in metres.}
#' \item{delta.elev}{delta elevation (metres).}
#' \item{VAM}{"vertical ascent metres per second".}
#' \item{power.W}{power readings (Watts).}
#' \item{power.smooth.W}{an exponentially-weighted 25-second moving average of power values.}
#' \item{work.kJ}{cumulative work (kilojoules).}
#' \item{Wexp.kJ}{W' expended in units of kilojoules. See \code{?Wbal} and references therein.}
#' \item{cadence.rpm}{pedalling cadence (revolutions per minute).}
#' \item{hr.bpm}{Heart rate (beats per minute).}
#' \item{lap}{a numeric vector of lap "levels". Will only have values > 1 if lap data is available.}
#' }
#' @name cycleRdata
NULL
