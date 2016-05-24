speedGPS <- function(xyz, time, rOut = c(.25, .75), naValue = .25) {
    speed <-  trackSpeed(xyz, time)
    rIn <- range(speed, na.rm = TRUE)
    value <- diff(rOut) * (speed - rIn[1])/diff(rIn) + rOut[1]
    ifelse(is.na(value), rep(naValue, length(value)), value)
}
