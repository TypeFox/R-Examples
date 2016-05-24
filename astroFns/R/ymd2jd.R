ymd2jd <-
function(yr=2012, mo=1, dy=1) {

# Year, month, date to Julian day conversion
# Rounds to nearest Julian date at input
# Input is normally three numerical values but the first argument may also
#    be a POSIXt class
# Reference time is 00:00 on date
#
# From Fliegel & Van Flandern, Comm. ACM 10, 657 (1968)
#   Note: algorithm uses FORTRAN integer mathematics
#   Also: Explanatory Supplement to the Astronomical Almanac
# A. Harris, U. Maryland Astronomy, 4/18/2008

    if (class(yr)[1]=='POSIXt') {
        # convert timezone to UTC, extract ymd
        yr <- as.POSIXlt(yr, tz='UTC')
        dy <- as.numeric(format(yr, '%d'))
        mo <- as.numeric(format(yr, '%m'))
        yr <- as.numeric(format(yr, '%Y'))
    }

    dy <- round(dy)
    mo <- round(mo)
    yr <- round(yr)

    jd <- dy - 32075 + trunc(1461*(yr + 4800 + trunc((mo - 14)/12))/4) +
          trunc(367*(mo - 2 - 12*trunc((mo-14)/12))/12) -
          trunc(3*(trunc((yr + 4900 + trunc((mo - 14)/12))/100))/4)

# convert from noon to 0h
jd - 0.5

}

