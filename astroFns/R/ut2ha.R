'ut2ha' <-
function(yr=2012, mo=1, dy=1, hr=0, mi=0, se=0,
         ra.sou='13h 31m 08.3s', lon.obs="W 79d 50m 23.4s") {

# Calculate hour angle from date, time, source RA, and observatory position.
# Default longitude for NRAO Green Bank.
#
# A. Harris 2012.05.18, 2012.07.05

    # times; test first for Sys.time() as input
    if (any(class(yr)=='POSIXct')) {
        # convert timezone to UTC, extract ymd
        yr <- as.POSIXlt(yr, tz='UTC')
        se <- as.numeric(format(yr, '%S'))
        mi <- as.numeric(format(yr, '%M'))
        hr <- as.numeric(format(yr, '%H'))
        dy <- as.numeric(format(yr, '%d'))
        mo <- as.numeric(format(yr, '%m'))
        yr <- as.numeric(format(yr, '%Y'))
    }

    # hour angle
    ha <- ut2lst(yr, mo, dy, hr, mi, se, lon.obs) - hms2rad(ra.sou)*12/pi
    # between -12h and 12h
    ha <- ha%%24
    idxp <- which(ha > 12)
    idxn <- which(ha < -12)
    ha[idxp] <- ha[idxp] - 24
    ha[idxn] <- ha[idxn] + 24

    ha  # inherits class fracHrs from ut2lst
}
