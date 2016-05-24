#' Distance between two locations
#' 
#' Distance between two locations set out by clicking on a geo--plot.
#' 
#' 
#' @param scale Unit of returned distance, default "nmi" for nautical miles,
#' all other values return kilometers
#' @param type Display points "p" or lines "l" between clicks, "n" for nothing
#' @return Returns distance between two point clicks on a geoplot.
#' @note Rather limited functionality, could be built further?
#' @seealso Calls \code{\link{arcdist}} and \code{\link{geolocator}}.
#' @keywords iplot
#' @export locdist
locdist <-
function(scale = "nmi", type = "p")
{
        lat <- geolocator(n = 2, type = type)
        x <- arcdist(lat$lat[1], lat$lon[1], lat$lat[2], lat$lon[2])
        if(scale == "km")
                x <- x * 1.852
        return(x)
}

