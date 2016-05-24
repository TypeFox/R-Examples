#' Add arrow to plot.
#' 
#' Adds arrow to a plot, given base and point position along with half angle.
#' 
#' 
#' @param pos List with components \code{lat, lon} of length 2 or more (only
#' first two used).
#' @param angle Angle (< 90) determining arrow sharpness
#' @param col Arrow color, defaults to color 2 of palette in effect.
#' @return Mostly used for plotting, but the function returns the arrow polygon
#' invisibly.
#' @seealso \code{\link{geopolygon}}, \code{\link{invProj}}
#' @keywords aplot
#' @examples
#' 
#' geoplot()
#' Arrow(list(lat=c(65,65.5),lon=c(-19, -18)),angle=60,col="brown")
#' Arrow(list(lat=c(65,65.5),lon=c(-19, -18)),angle=45,col="red")
#' Arrow(list(lat=c(65,65.5),lon=c(-19, -18)),angle=30,col="green")
#' Arrow(list(lat=c(65,65.5),lon=c(-19, -18)),angle=30,col="blue")
#' 
#' @export Arrow
Arrow <-
function (pos, angle = 15, col = 2)
{
    pos <- Proj(pos)
    dx <- -diff(pos$x)
    dy <- -diff(pos$y)
    d <- sqrt(dy * dy + dx * dx)
    d1 <- d * tan((angle * pi)/180)
    p1y <- pos$y[1] + d1/d * dx
    p1x <- pos$x[1] - d1/d * dy
    p2y <- pos$y[1] - d1/d * dx
    p2x <- pos$x[1] + d1/d * dy
    d <- data.frame(y = c(pos$y[2], p1y, p2y, pos$y[2]), x = c(pos$x[2],
        p1x, p2x, pos$x[2]))
    d <- invProj(d)
    geopolygon(d,col=col)
    return(invisible(data.frame(lat = d$lat, lon = d$lon)))
}

