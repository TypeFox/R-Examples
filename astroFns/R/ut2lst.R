ut2lst <- function(yr=2012, mo=1, dy=1, hr=0, mi=0, se=0,
                   lon.obs="W 79d 50.5m") {

# Calculate local sidereal time from date, time, and observatory position
# Ref (especially for the equation of the equinoxes):
#    http://aa.usno.navy.mil/faq/docs/GAST.php
# but with more precise calculation of GMST1.
# Spot checks show values match tabulated values in The Astronomical Almanac
# within ~0.01 seconds.
# Longitude for NRAO Green Bank.
#
# A. Harris 2012.05.17, 2012.07.05

    # longitude in hours
    lon.h <- dms2rad(lon.obs)*12/pi

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

    h <- (hr + (mi + (se/60))/60)
    d <- ymd2jd(yr, mo, dy) + h/24 - 2451545.0
    eqeq <- (-0.000319*sin((125.04 - 0.052954*d)*pi/180) -
               0.000024*sin((280.47 + 0.98565*d)*pi/90)) *
                   cos((23.4393 - 0.0000004*d)*pi/180)
    h <- 1.002737909350795 * h  # first-order (constant day length) term only

    lst <- gmst1(yr, mo, dy) + h + lon.h + eqeq
    lst <- (lst + 24)%%24  # ensure always between 0 and 24h

    class(lst) <- 'fracHrs'
    lst
}
