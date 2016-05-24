
setClassUnion("DateTime", c("POSIXt", "POSIXct", "POSIXlt"))

setClass("GPSTrack",
  representation(latitude = "numeric", longitude = "numeric",
                 elevation = "numeric", time = "DateTime")
)

scanGPSTrack <- function(con, 
    fields = list(date = "", time = "",
                    lat = 0., lon = 0., el = 0.),
    dateTimeFormat = "%Y-%m-%d %H:%M:%S") {
    data <- scan(con, fields)
    txt <- textConnection(paste(data$date, data$time), "r")
    dateTime <- .scanDateTime(txt, dateTimeFormat)
    close(txt)
    new("GPSTrack", latitude = data$lat, longitude = data$lon,
          elevation = data$el, time = dateTime)
}

.scanDateTime <- function(con, dateTimeFormat) {
    as.POSIXct(strptime(readLines(con), dateTimeFormat))
}
