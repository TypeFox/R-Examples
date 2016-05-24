`gmst1` <-
function(yr=2012, mo=1, dy=1) {
# Greenwich mean sidereal time (GMST) at 0h UT1
# Ref: Explanatory Supplement to the Astronomical Almanac
# Seidelmann (ed), c. 1992
# also Aoki, S. et al. 2008, Astr. Ap. 105, 359-361
#
# Requires Julian date conversion routine ymd2jd
#
# Output is in fractional day [0-1] with class fracHrs
#
# Convert date to JD at 0h (0h is JD - 0.5)
# For fractional day multiply UT1 fractional day by
# 1.002737909350795 to get fractional sidereal day,
# ignoring small higher order corrections (2.241-4).
# Mean sidereal day is 23h56m04.090542 of UT1
#
# AH 2010.02.11-2012.06.15

    if (class(yr)[1]=='POSIXt') {
        # convert timezone to UTC, extract ymd
        yr <- as.POSIXlt(yr, tz='UTC')
        dy <- as.numeric(format(yr, '%d'))
        mo <- as.numeric(format(yr, '%m'))
        yr <- as.numeric(format(yr, '%Y'))
    }

    # Calculation
    du <- ymd2jd(yr, mo, dy) - 2451545
    Tu <-  du/36525
    # GMST at 0h in seconds (2.24-1)
    gmst1 <- 24110.54841 + (8640184.812866 + (0.093104 - 6.2e-6*Tu)*Tu)*Tu
    # In hours
    gmst1 <- (gmst1/86400)%%1 * 24

    class(gmst1) <- 'fracHrs'
    gmst1
}

