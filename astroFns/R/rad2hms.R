rad2hms <-
function(rad=1, places=1) {
    # Check places for rounding
    places <- round(places,  0)
    if (places <= 0) places <- 0
    if (places > 6) places <- 6
    # Convert radians to hms
    sgn <- ifelse(rad>=0, 1, -1)
    hrs <- rad*12/pi
    ahrs <- abs(hrs%%24)
    # Work out whole hours, minutes, fractional seconds
    hrsfrac <- ahrs%%1
    min <- (hrsfrac%%1)*60
    sec <- round((min%%1)*60, places)
    # Clean up rounding
    for (i in length(sec)) {
        if(sec[i]>=60) {
            sec[i] <- sec[i]-60
            min[i] <- min[i]+1
        }
        if(min[i]>=60) {
            min[i] <- min[i]-60
            ahrs[i] <- ahrs[i]+1
        }
    }
    # Format output
    if (places == 0) {
        fmt <- '%02d:%02d:%02d'
    } else {
        pl <- 3+places + (places)/10
        fmt <- paste('%02d:%02d:%0#', pl, 'f', sep='')
    }
    sprintf(fmt, trunc(ahrs), trunc(min), sec)
}

