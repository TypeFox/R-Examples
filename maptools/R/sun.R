".deg" <- function(radian) 180 * radian / pi
".rad" <- function(degree) pi * degree / 180

".julianD" <- function(year, month, day)
{
    ## Value: Numeric Julian day without fractions.
    ## --------------------------------------------------------------------
    ## Arguments: year=4-digit year, month=1-12, day=1-31, all integers.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    corr <- month <= 2
    year[corr] <- year[corr] - 1
    month[corr] <- month[corr] + 12
    a <- floor(year / 100)
    b <- 2 - a + floor(a / 4)
    floor(365.25 * (year + 4716)) + floor(30.6001 * (month + 1)) +
        day + b - 1524.5
}

## Takes Julian day and returns centuries since J2000.0
".cent2000JD" <- function(jd) (jd - 2451545) / 36525

## Takes number of centuries since J2000.0 and returns julian day
".julianD2000" <- function(jc) jc * 36525 + 2451545

".geomMeanLonSun" <- function(jc)
{
    ## Value: The geometric mean longitude of the sun in degrees.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    l0 <- 280.46646 + jc * (36000.76983 + 0.0003032 * jc)
    updown <- !is.finite(l0) | l0 > 360 | l0 < 0
    l0[updown] <- l0[updown] %% 360
    l0
}

".geomMeanAnomSun" <- function(jc)
{
    ## Value: Numeric, geometric mean anomaly of the sun in degrees.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    357.52911 + jc * (35999.05029 - 1.537e-4 * jc)
}

".eccentrEarthOrb" <- function(jc)
{
    ## Value: Numeric, unitless eccentricity of the Earth's orbit.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    0.016708634 - jc * (4.2037e-5 + 1.267e-7 * jc)
}

".sunEqCenter" <- function(jc)
{
    ## Value: Numeric, position of the center of the sun in degrees.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    m <- .geomMeanAnomSun(jc)
    mrad <- .rad(m)
    sinm <- sin(mrad)
    sin2m <- sin(mrad * 2)
    sin3m <- sin(mrad * 3)
    sinm * (1.914602 - jc * (4.817e-3 + 1.4e-5 * jc)) + sin2m *
        (1.9993e-2 - 1.01e-4 * jc) + sin3m * 2.89e-4
}

".sunTrueLon" <- function(jc)
{
    ## Value: Numeric, sun's true longitude in degrees.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    l0 <- .geomMeanLonSun(jc)
    eqc <- .sunEqCenter(jc)
    l0 + eqc
}

".sunTrueAnom" <- function(jc)
{
    ## Value: Numeric, sun's true anomaly in degrees.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    m <- .geomMeanAnomSun(jc)
    eqc <- .sunEqCenter(jc)
    m + eqc
}

".sunRadVec" <- function(jc)
{
    ## Value: Numeric, sun's radius vector in AUs.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    v <- .sunTrueAnom(jc)
    eo <- .eccentrEarthOrb(jc)
    (1.000001018 * (1 - eo * eo)) / (1 + eo * cos(.rad(v)))
}

".sunApparentLon" <- function(jc)
{
    ## Value: Numeric, sun's apparent longitude in degrees
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    tl <- .sunTrueLon(jc)
    om <- 125.04 - 1934.136 * jc
    tl - 0.00569 - 0.00478 * sin(.rad(om))
}

".meanObliqEcliptic" <- function(jc)
{
    ## Value: Numeric, mean obliquity of the ecliptic in degrees.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    s <- 21.448 - jc * (46.815 + jc * (5.9e-4 - jc * 0.001813))
    23 + (26 + (s / 60)) / 60
}

".obliqCorr" <- function(jc)
{
    ## Value: Numeric, the corrected obliquity of the ecliptic in degrees.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    e0 <- .meanObliqEcliptic(jc)
    om <- 125.04 - 1934.136 * jc
    e0 + 0.00256 * cos(.rad(om))
}

".sunRtAscension" <- function(jc)
{
    ## Value: Numeric, the sun's right ascension in degrees.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    e0ok <- .obliqCorr(jc)
    la <- .sunApparentLon(jc)
    tananum <- cos(.rad(e0ok)) * sin(.rad(la))
    tanadenom <- cos(.rad(la))
    .deg(atan2(tananum, tanadenom))
}

".sunDeclination" <- function(jc)
{
    ## Value: Numeric, sun's declination in degrees.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    e0ok <- .obliqCorr(jc)
    la <- .sunApparentLon(jc)
    sinjc <- sin(.rad(e0ok)) * sin(.rad(la))
    .deg(asin(sinjc))
}

".eqTime" <- function(jc)
{
    ## Value: Numeric, equation of the difference between true solar and
    ## mean solar times.
    ## --------------------------------------------------------------------
    ## Arguments: jc=number of centuries since J2000.0
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    epsi <- .obliqCorr(jc)
    l0 <- .geomMeanLonSun(jc)
    ecc <- .eccentrEarthOrb(jc)
    m <- .geomMeanAnomSun(jc)
    y <- tan(.rad(epsi) / 2) ^ 2
    sin2l0 <- sin(2 * .rad(l0))
    sinm <- sin(.rad(m))
    cos2l0 <- cos(2 * .rad(l0))
    sin4l0 <- sin(4 * .rad(l0))
    sin2m <- sin(2 * .rad(m))
    etime <- y * sin2l0 - 2 * ecc * sinm + 4 * ecc * y * sinm * cos2l0 -
        0.5 * y * y * sin4l0 - 1.25 * ecc * ecc * sin2m
    .deg(etime) * 4
}

".hangleCrepuscule" <- function(lat, solarDec, solarDep,
                                direction=c("dawn", "dusk"))
{
    ## Value: Numeric, hour angle of the sun at dawn or dusk in radians.
    ## --------------------------------------------------------------------
    ## Arguments: solarDec=declination angle of the sun in degrees;
    ## solarDep=angle of the sun below the horizon in degrees;
    ## dawn=logical indicating whether dawn or dusk hour angle should be
    ## returned.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    latrad <- .rad(lat)
    sdrad <- .rad(solarDec)
    haarg <- (cos(.rad(90 + solarDep)) / (cos(latrad) * cos(sdrad)) -
              tan(latrad) * tan(sdrad))
    haarg[abs(haarg) >= 1] <- NA
    angle <- acos(haarg)
    switch(direction, dawn=angle, dusk=-angle)
}

".hangleSunriset" <- function(lat, solarDec, direction=c("sunrise", "sunset"))
{
    ## Value: Numeric, hour angle of the sun at sunrise or sunset in
    ## radians.
    ## --------------------------------------------------------------------
    ## Arguments: lat=numeric, latitude of observer in degrees;
    ## solarDec=declination angle of the sun in degrees;
    ## sunrise=logical indicating whether sunrise or sunset hour angle
    ## should be returned.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    latrad <- .rad(lat)
    sdrad <- .rad(solarDec)
    haarg <- (cos(.rad(90.833)) / (cos(latrad) * cos(sdrad)) -
              tan(latrad) * tan(sdrad))
    haarg[abs(haarg) >= 1] <- NA
    angle <- acos(haarg)
    switch(direction, sunrise= angle, sunset=-angle)
}

".crepusculeUTC" <- function(jd, lon, lat, solarDep,
                             direction=c("dawn", "dusk"))
{
    ## Value: Numeric, UTC time of dawn or dusk, in minutes from zero Z.
    ## --------------------------------------------------------------------
    ## Arguments: jd=julian day (real);
    ## lon=lat=longitude and latitude, respectively, of the observer in
    ## degrees;
    ## solarDep=angle of the sun below the horizon in degrees;
    ## dawn=logical indicating whether dawn or dusk UTC should be
    ## returned.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    jc <- .cent2000JD(jd)
    eqtime <- .eqTime(jc)
    solarDec <- .sunDeclination(jc)
    switch(direction,
           dawn = hourangle <- .hangleSunriset(lat, solarDec, direction="sunrise"),
           dusk = hourangle <- .hangleSunriset(lat, solarDec, direction="sunset"))
    d <- lon - .deg(hourangle)
    tdiff <- 4 * d
    timeUTC <- 720 + tdiff - eqtime
    newt <- .cent2000JD(.julianD2000(jc) + timeUTC / 1440)
    eqtime <- .eqTime(newt)
    solarDec <- .sunDeclination(newt)
    switch(direction,
           dawn = {
               hourangle <- .hangleCrepuscule(lat, solarDec, solarDep,
                                              direction="dawn")},
           dusk = {
               hourangle <- .hangleCrepuscule(lat, solarDec, solarDep,
                                              direction="dusk")})
    d <- lon - .deg(hourangle)
    tdiff <- 4 * d
    720 + tdiff - eqtime
}

".sunrisetUTC" <- function(jd, lon, lat, direction=c("sunrise", "sunset"))
{
    ## Value: Numeric, UTC time of sunrise or sunset, in minutes from zero
    ## Z.
    ## --------------------------------------------------------------------
    ## Arguments: jd=julian day (real);
    ## lon=lat=longitude and latitude, respectively, of the observer in
    ## degrees;
    ## sunrise=logical indicating whether sunrise or sunset UTC should be
    ## returned.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    jc <- .cent2000JD(jd)
    eqtime <- .eqTime(jc)
    solarDec <- .sunDeclination(jc)
    switch(direction,
           sunrise = hourangle <- .hangleSunriset(lat, solarDec, direction="sunrise"),
           sunset = hourangle <- .hangleSunriset(lat, solarDec, direction="sunset"))
    d <- lon - .deg(hourangle)
    tdiff <- 4 * d
    timeUTC <- 720 + tdiff - eqtime
    newt <- .cent2000JD(.julianD2000(jc) + timeUTC / 1440)
    eqtime <- .eqTime(newt)
    solarDec <- .sunDeclination(newt)
    switch(direction,
           sunrise = hourangle <- .hangleSunriset(lat, solarDec, direction="sunrise"),
           sunset = hourangle <- .hangleSunriset(lat, solarDec, direction="sunset"))
    d <- lon - .deg(hourangle)
    tdiff <- 4 * d
    720 + tdiff - eqtime
}

".redoLonLat" <- function(lon, lat)
{
    ## Value: Matrix of latitude and longitude with + N, + W values, and
    ## with latitudes < -89.8 or > 89.8 fixed to -89.8 and 89.8,
    ## respectively.
    ## --------------------------------------------------------------------
    ## Arguments: lon=lat=longitude and latitude, respectively, in degrees
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    newlon <- -lon
    newlat <- lat
    newlat[newlat > 89.8] <- 89.8
    newlat[newlat < -89.8] <- -89.8
    cbind(newlon, newlat)
}

".timeData" <- function(time)
{
    ## Value: list with numeric vectors year, month, day, offset hours
    ## from GMT, and whether day light savings is in effect.
    ## --------------------------------------------------------------------
    ## Arguments: time=POSIXct
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    time.gmt <- as.POSIXct(format(time), tz="GMT")
    time.plt <- as.POSIXlt(time)
    timezone <- as.numeric(difftime(time.gmt, time, units="hours"))
    year <- as.integer(format(time.plt, "%Y"))
    month <- as.integer(format(time.plt, "%m"))
    day <- as.integer(format(time.plt, "%d"))
    hour <- as.integer(format(time.plt, "%H"))
    min <- as.integer(format(time.plt, "%M"))
    sec <- as.integer(format(time.plt, "%S"))
    list(year=year, month=month, day=day, hour=hour, min=min, sec=sec,
         timezone=timezone, dlstime=0, tz=attr(time, "tzone"))
}

".crepuscule" <- function(lon, lat, year, month, day, timezone,
                          dlstime, solarDep, direction=c("dawn", "dusk"))
{
    ## Value: Numeric, time of dawn in local time (days)
    ## --------------------------------------------------------------------
    ## Arguments: lon=lat=longitude and latitude of the observer in degrees;
    ## year=4-digit year; month=1-12; day=1-31;
    ## timezone=time zone hours shift relative to UTC (hours);
    ## dlstime=1 or 0 to indicate daylight savings time or not;
    ## solarDep=angle of the sun below the horizon in degrees;
    ## dawn=logical to indicate whether dawn or dusk time should be
    ## returned
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    ll <- .redoLonLat(lon, lat)
    jd <- .julianD(year, month, day)
    switch(direction,
           dawn = {
               risetTimeGMT <- .crepusculeUTC(jd, ll[, 1], ll[, 2], solarDep,
                                              direction="dawn")},
           dusk = {
               risetTimeGMT <- .crepusculeUTC(jd, ll[, 1], ll[, 2], solarDep,
                                              direction="dusk")})
    risetTimeLST <- risetTimeGMT + (60 * timezone) + (dlstime * 60)
    risetTimeLST / 1440
}

".sunriset" <- function(lon, lat, year, month, day, timezone,
                        dlstime, direction=c("sunrise", "sunset"))
{
    ## Value: Numeric, time of sunrise in local time (days)
    ## --------------------------------------------------------------------
    ## Arguments: lon=lat=longitude and latitude of the observer in degrees;
    ## year=4-digit year; month=1-12; day=1-31;
    ## timezone=time zone hours shift relative to UTC (hours);
    ## dlstime=1 or 0 to indicate daylight savings time or not;
    ## sunrise=logical to indicate whether sunrise or sunset time should
    ## be returned
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    ll <- .redoLonLat(lon, lat)
    jd <- .julianD(year, month, day)
    switch(direction,
           sunrise = {
               risetTimeGMT <- .sunrisetUTC(jd, ll[, 1], ll[, 2],
                                            direction="sunrise")},
           sunset = {
               risetTimeGMT <- .sunrisetUTC(jd, ll[, 1], ll[, 2],
                                            direction="sunset")})
    risetTimeLST <- risetTimeGMT + (60 * timezone) + (dlstime * 60)
    risetTimeLST / 1440
}

".solarnoon" <- function(lon, lat, year, month, day, timezone, dlstime)
{
    ## Value: Numeric, time of solar noon in local time (days)
    ## --------------------------------------------------------------------
    ## Arguments: lon=lat=longitude and latitude of the observer in degrees;
    ## year=4-digit year; month=1-12; day=1-31;
    ## timezone=time zone hours shift relative to UTC (hours);
    ## dlstime=1 or 0 to indicate daylight savings time or not;
    ## solarDep=angle of the sun below the horizon in degrees
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    ll <- .redoLonLat(lon, lat)
    jd <- .julianD(year, month, day)
    jc <- .cent2000JD(jd)
    newt <- .cent2000JD(.julianD2000(jc) + 0.5 + ll[, 1] / 360)
    eqtime <- .eqTime(newt)
    ## solarNoonDec <- .sunDeclination(newt)
    solNoonUTC <- 720 + (ll[, 1] * 4) - eqtime
    solarnoon <- solNoonUTC + (60 * timezone) + (dlstime * 60)
    solarnoon / 1440
}

".solarpos" <- function(lon, lat, year, month, day, hours, minutes,
                        seconds, timezone, dlstime)
{
    ## Value: matrix with solar azimuth (in degrees from N) and solar
    ##  elevation.
    ##  --------------------------------------------------------------------
    ##  Arguments: lon=lat=longitude and latitude of the observer in
    ##  degrees; year=4-digit year; month=1-12; day=1-31; hours=0-23;
    ##  minutes=0-59; seconds=0-59; timezone=time zone hours shift
    ##  relative to UTC (hours); dlstime=1 or 0 to indicate daylight
    ##  savings time or not;
    ##  --------------------------------------------------------------------
    ##  Author: Sebastian Luque
    ##  --------------------------------------------------------------------
    ll <- .redoLonLat(lon, lat)
    zone <- -timezone
    hh <- hours - ((dlstime * 60) / 60)
    timenow <- hh + minutes / 60 + seconds / 3600 + zone
    jd <- .julianD(year, month, day)
    jc <- .cent2000JD(jd + timenow / 24)
    ## earthRadVec <- r <- .sunRadVec(jc)
    solarDec <- .sunDeclination(jc)
    eqtime <- .eqTime(jc)
    solarTimefix <- eqtime - 4 * ll[, 1] + 60 * zone
    trueSolarTime <- hh * 60 + minutes + seconds / 60 + solarTimefix
    corrsol <- trueSolarTime > 1440
    trueSolarTime[corrsol] <- trueSolarTime[corrsol] %% 1440
    hourangle <- trueSolarTime / 4 - 180
    hourangle[hourangle < -180] <- hourangle[hourangle < -180] + 360
    harad <- .rad(hourangle)
    csz <- sin(.rad(ll[, 2])) * sin(.rad(solarDec)) + cos(.rad(ll[, 2])) *
        cos(.rad(solarDec)) * cos(harad)
    csz[csz > 1] <- 1
    csz[csz < -1] <- -1
    zenith <- .deg(acos(csz))
    azDenom <- cos(.rad(ll[, 2])) * sin(.rad(zenith))
    azimuth <- numeric(length(azDenom))
    hiaD <- abs(azDenom) > 0.001        # if (Abs(azDenom) > 0.001) ... BEG
    azRad <- ((sin(.rad(ll[hiaD, 2])) * cos(.rad(zenith[hiaD]))) -
              sin(.rad(solarDec))) / azDenom[hiaD]
    zz <- abs(azRad) > 1                # if (Abs(azRad) > 1) ... BEG
    azRad[zz & azRad < 0] <- -1         # if (azRad < 0) ... BEG
    azRad[zz & !azRad < 0] <- 1         # if (azRad < 0) ... END
                                        # if (Abs(azRad) > 1) ... END
    azimuth1 <- 180 - .deg(acos(azRad))
    zz <- hourangle[hiaD] > 0           # if (hourangle > 0) ... BEG
    azimuth1[zz] <- -azimuth1[zz]       # if (hourangle > 0) ... END
    azimuth[hiaD] <- azimuth1
    loaD <- !hiaD
    azimuth[loaD & ll[, 2] > 0] <- 180  # if (latitude > 0) ... BEG
    azimuth[loaD & !ll[, 2] > 0] <- 0   # if (latitude > 0) ... END
                                        # if (Abs(azDenom) > 0.001) ... END
    azimuth[azimuth < 0] <- azimuth[azimuth < 0] + 360
    exoatmEl <- 90 - zenith
    refracCorr <- numeric(length(exoatmEl))
    hiR <- exoatmEl > 85                # if (exoatmElevation > 85) ... BEG
    refracCorr[hiR] <- 0
    loR <- !hiR
    zz <- loR & exoatmEl > 5
    te <- tan(.rad(exoatmEl[zz]))
    refracCorr[zz] <- 58.1 / te - 0.07 / (te^3) + 8.6e-5 / te^5
    zz <- loR & !exoatmEl > 5 & exoatmEl > -0.575
    step1 <- -12.79 + exoatmEl[zz] * 0.711
    step2 <- 103.4 + exoatmEl[zz] * step1
    step3 <- -518.2 + exoatmEl[zz] * step2
    refracCorr[zz] <- 1735 + exoatmEl[zz] * step3
    zz <- loR & !exoatmEl > 5 & !exoatmEl > -0.575
    te <- tan(.rad(exoatmEl[zz]))
    refracCorr[zz] <- -20.774 / te
    refracCorr <- refracCorr / 3600
    solarzen <- zenith - refracCorr
    cbind(azimuth=azimuth, elevation=90 - solarzen)
}
