if(!isClass("GPSTrack")) {
  source("./DateTimeTools.R")
  source("./GPSTrackClass.R")
}

scanGPSTrack <- function(con,  ...) {
    fields <- list(date = "", time = "",
                    x = 0., y = 0., z = 0.)
    data <- scan(con, fields)
    dateTime <- as.POSIXct(strptime(paste(date, time),
                                    "20%y-%m-%d %H:%M:%S"))
    new("GPSTrack", latitude = data$x, longtitude = data$y,
          elevation = data$z, time = dateTime)
}
