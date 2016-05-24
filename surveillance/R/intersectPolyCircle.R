################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Compute the intersection of a circular domain with a polygonal domain of
### various classes (currently: owin, gpc.poly, or SpatialPolygons)
###
### Copyright (C) 2009-2015 Sebastian Meyer
### $Revision: 1169 $
### $Date: 2015-01-13 16:05:56 +0100 (Die, 13. Jan 2015) $
################################################################################


intersectPolyCircle.gpc.poly <- function (object, center, radius,
                                          npoly = 32, useGEOS = FALSE, ...)
{
    if (useGEOS) {
        loadNamespace("rgeos")  # coerce gpc.poly to SpatialPolygons
        res <- intersectPolyCircle.SpatialPolygons(as(object, "SpatialPolygons"),
                                               center, radius, npoly)
        as(res, "gpc.poly")  # also defined in rgeos
    } else {
        gpclibCheck()
        circle <- discpoly(center, radius, npoly = npoly, class = "gpc.poly")
        gpclib::intersect(circle, object)  # this order seems to be faster
    }
}

intersectPolyCircle.owin <- function (object, center, radius, npoly = 32, ...)
{
    circle <- disc(radius = radius, centre = center, npoly = npoly)
    res <- intersect.owin(circle, object)  # order does not affect runtime
    ## ensure "polygonal" type (because of rescue.rectangle in intersect.owin)
    as.polygonal(res)
}

intersectPolyCircle.SpatialPolygons <- function (object, center, radius,
                                                 npoly = 32, ...)
{
    circle <- discpoly(center, radius, npoly = npoly, class = "Polygon")
    circleSpP <- SpatialPolygons(list(Polygons(list(circle), "0")))
    ## ensure that circleSpP has exactly the same proj4string as 'object'
    circleSpP@proj4string <- object@proj4string
    rgeos::gIntersection(circleSpP, object)
}
