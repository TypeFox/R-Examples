marinus <- function(lat, lon) {

################################################################################
# Function: marinus
# Purpose: Convert coordinates from latitude/longitude to the equidistant,
#          cylindric map projection
# Programmer: Denis White
# Date: September 5, 2002
# Last Revised: April 4, 2007
# Description:
#   This function converts x,y coordinates measured in units of latitude and 
#   longitude, i.e., geographic coordinates measured in decimal degrees, to 
#   coordinates in the equidistant, cylindric map projection measured in units of
#   kilometers.  The projection center is defined as the midpoint in latitude-
#   longitude space.  The map projection is here named after Marinus of Tyre
#   (see  J.P. Snyder. USGS Prof Paper 1395, p. 90).
# Arguments:
#   lat = vector of latitudes.
#   lon = vector of longitudes.
# Results:
#   A matrix with column names "x" and "y" containing the x and y coordinates 
#   in the equidistant, cylindric map projection measured in units of 
#   kilometers.
# Other Functions Required:  None
# Examples:
#   lat <- 45 + runif(100, -5, 5)
#   lon <- 120 + runif(100, -10, 10)
#   marinus(lat, lon)
################################################################################

    R <- 6371 # authalic radius of Clarke 1866 rounded to km
    rlat <- range(lat, na.rm=TRUE)
    rlon <- range(lon, na.rm=TRUE)
    clat <- diff(rlat) / 2 + rlat[1]
    clon <- diff(rlon) / 2 + rlon[1]
    x <- R * (lon-clon) * pi/180 * cos(clat*pi/180)
    y <- R * (lat-clat) * pi/180
    x[is.na(lat)] <- lon[is.na(lat)]
    xy <- data.frame(x=x, y=y)

    xy
}
