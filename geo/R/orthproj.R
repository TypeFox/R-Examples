#' Performs an orthogonal projection to a curve.
#' 
#' Finds curve coordinates of points, given points and a curve it returns the
#' length of the points orhogonal projection and the distance traveled to the
#' projection from a given origin on the curve.
#' 
#' 
#' @param pts A list containing the points as \$lat and \$lon, you may also use
#' geolocator.
#' @param curve The curve to be used for projection.
#' @return Returns a dataframe with vectors pardist (length of orthogonal
#' projection) and perdist (how far traveled alongst the curve).
#' @section Side Effects: None.
#' @seealso \code{\link{geocurve}}, \code{\link{geolocator}}.
#' @examples
#' 
#'        
#' \dontrun{       geoplot(my.curve)     # Plot curve and initialize plot.
#'        geocurve(geolocator(),my.curve)         # Mark points.
#' }
#' @export orthproj
orthproj <-
function(pts, curve)
{
        pts1 <- lambert(pts$lat, pts$lon, 65, -18, 65)
        curve1 <- lambert(curve$lat, curve$lon, 65, -18, 65)
        pts$x <- pts1$x/1.852
        pts$y <- pts1$y/1.852
        curve$x <- curve1$x/1.852
        curve$y <- curve1$y/1.852
        pardist <- perdist <- rep(0, length(pts$lat))
        x <- .C("Curvedist", PACKAGE = "geo", 
                as.double(curve$x),
                as.double(curve$y),
                as.double(curve$dist),
                as.integer(length(curve$x)),
                as.double(pts$x),
                as.double(pts$y),
                as.double(perdist),
                as.double(pardist),
                as.integer(length(pts$y)))
        pardist <- x[[7]]
        perdist <- x[[8]]
        # Points inside curve get negative perdist.
        i <- geoinside(pts, reg = curve, option = 0, robust = F)
        perdist[i] <-  - perdist[i]
        return(list(pardist = pardist, perdist = perdist))
}

