
setOldClass("POSIXt") ## a lie because POSIXt really has an illegal inheritance
setClassUnion("DateTime", "POSIXt")

setClass("GPSTrack",
  representation(latitude = "numeric", longitude = "numeric",
                 elevation = "numeric", time = "DateTime")
)

scanGPSTrack <- function(con,  ...) {
    fields <- list(date = "", time = "",
                    x = 0., y = 0., z = 0.)
    data <- scan(con, fields)
    dateTime <- scanDateTime(textConnection(
                  paste(data$date, data$time)))
    new("GPSTrack", latitude = data$x, longitude = data$y,
          elevation = data$z, time = dateTime)
}

## simple function for the moment, suitable to read ISO time/date
## should eventually intuit time format via Perl?

scanDateTime <- function(con) {
    as.POSIXct(strptime(readLines(con), "20%y-%m-%d %H:%M:%S"))
}

geoCoords <- function(latitude, longitude, origin = c(40.7, -74.4)) {
  n = length(latitude)
  if(length(longitude) != n)
    stop("required equal length of latitude, longitude, got ", n,
  ", ", length(longitude))
  x = geoDist(latitude, longitude, latitude, rep(origin[[2]], n))
  *sign(longitude-origin[[2]])
  y = geoDist(latitude, longitude, rep(origin[[1]], n), longitude) *
  sign(latitude - origin[[1]])
  list(x=x, y=y)
}
