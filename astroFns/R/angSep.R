angSep <-
function(ra1, dec1, ra2, dec2) {

# Function calculates angular separation of two sources
# Input in RA and Dec, hms and dms, comma-separated strings
# Returns angular distance in decimal degrees
#
# A. Harris 2010.2.14

    # Convert RA, Dec to radians
    ra1 <- hms2rad(ra1)
    dec1 <- dms2rad(dec1)
    ra2 <- hms2rad(ra2)
    dec2 <- dms2rad(dec2)
    # Compute angular distance between two positions and return
    a <- acos(sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra1-ra2))
    a*180/pi

}

