station <- function(lon = NULL, lat = NULL, pch = 24, col = "blue",
    bg = "red", ...)
{    
    #
    # draw a world map
    #
    map.world(main = "World Map with the Distribution of the Stations")

    if (missing(lon) || missing(lat))
    {
        ## read in the station longitude-latitude data and do a conversion.
        lon <- WTdiff$lon*pi/180; lat <- WTdiff$lat*pi/180
    }
    else if (length(lon) != length(lat))
      stop("The longitudes and latitudes vector must have the same length!\n")
    else
    {
        lon <- lon*pi/180; lat <- lat*pi/180
    }
       
    #
    # add the symbols for all the stations.
    #
    points(lon, lat, pch = pch, col = col, bg = bg, ...)

}


