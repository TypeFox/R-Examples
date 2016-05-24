`elev` <-
function(dec.sou = "33d 09m 35.0s", ha = 0, lat.obs = "38d 25m 59.2s") {
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
    # Calculate elevation in degrees
    asin(sin(lat.obs)*sin(dec.sou)+cos(lat.obs)*cos(dec.sou)*cos(ha))*180/pi
}

