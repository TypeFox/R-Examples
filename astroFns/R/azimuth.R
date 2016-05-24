`azimuth` <-
function(dec.sou='33d 09m 35.0s', ha=0, lat.obs='38d 25m 59.2s') {
    # Function calculates elevation [deg] given source dec.sou, HA
    # d is source dec.sou [dms string]
    # h is hour angle [dec.souimal hours]
    # lat.obs is observatory latitude [dms string]
    # "Astrophysical Formulae," K.R. Lang, c. 1986, 5-45
    # AH 12/28/08, 2/22/10

    # Convert observatory latitude, source dec.sou and hour angle to radians
    dec.sou <- dms2rad(dec.sou)
    ha <- unclass(ha)*pi/12
    lat.obs <- dms2rad(lat.obs)
    # Calculate elevation and parts for azimuth
    a <- asin(sin(lat.obs)*sin(dec.sou)+cos(lat.obs)*cos(dec.sou)*cos(ha))
    az1 <- -cos(dec.sou)*sin(ha)/cos(a)
    az2 <- (sin(dec.sou)*cos(lat.obs)-cos(lat.obs)*cos(dec.sou)*cos(ha))/cos(a)
    # Then azimuth
    az <- atan2(az1, az2) * 180/pi
    idx <- which(az<0)
    az[idx] <- az[idx]+360
    az
}

