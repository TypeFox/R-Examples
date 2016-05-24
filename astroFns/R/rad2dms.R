rad2dms <-
function(rad=1, places=2) {
    # Check places for rounding
    places <- round(places,  0)
    if (places <= 0) places <- 0
    if (places > 6) places <- 6
    # Convert radians to dms
    sgn <- ifelse(rad>=0, 1, -1)
    deg <- rad*180/pi
    adeg <- abs(deg)
    # Work out whole degrees, minutes, fractional seconds
    degfrac <- adeg%%1
    min <- (degfrac%%1)*60
    sec <- round((min%%1)*60, places)
    # Clean up rounding
    for (i in length(sec)) {
        if(sec[i]>=60) {
            sec[i] <- sec[i]-60
            min[i] <- min[i]+1
        }
        if(min[i]>=60) {
            min[i] <- min[i]-60
            adeg[i] <- adeg[i]+1
        }
    }
    # Format output
    if (places == 0) {
        fmt <- '%s%1d:%02d:%02d'
    } else {
        pl <- 3+places + (places)/10
        fmt <- paste('%s%1d:%02d:%0#', pl, 'f', sep='')
    }
    sprintf(fmt, ifelse(sgn<0, '-', '+'), trunc(adeg), trunc(min), sec)
}


