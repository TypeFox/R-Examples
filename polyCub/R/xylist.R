################################################################################
### Part of the R package "polyCub".
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Copyright (C) 2012-2014 Sebastian Meyer
### Time-stamp: <[xylist.R] 2015-02-25 21:18 (CET) by SM>
################################################################################


##' Convert Various Polygon Classes to a Simple List of Vertices
##' 
##' Different packages concerned with spatial data use different polygon
##' specifications, which sometimes becomes very confusing (see Details below).
##' To be compatible with the various polygon classes, package \pkg{polyCub}
##' uses an S3 class \code{"xylist"}, which represents
##' polygons by their core feature only, a list of lists of vertex coordinates
##' (see the "Value" section below).
##' The generic function \code{xylist} can deal with the
##' following polygon classes:
##' \itemize{
##' \item{\code{"\link[=owin.object]{owin}"} from package \pkg{spatstat}}
##' \item{\code{"\link[rgeos:gpc.poly-class]{gpc.poly}"} from package
##' \pkg{rgeos} (or \pkg{gpclib})}
##' \item{\code{"\linkS4class{Polygons}"} from package \pkg{sp}
##' (as well as \code{"\linkS4class{Polygon}"} and
##' \code{"\linkS4class{SpatialPolygons}"})}
##' }
##' The (somehow useless) default \code{xylist}-method
##' does not perform any transformation but only ensures that the polygons are
##' not closed (first vertex not repeated).
##' 
##' Different packages concerned with spatial data use different polygon
##' specifications with respect to:
##' \itemize{
##' \item{do we repeat the first vertex?}
##' \item{which direction represents holes?}
##' }
##' Package overview:
##' \describe{
##' \item{\pkg{sp}:}{\emph{Repeat} first vertex at the end (closed),
##' anticlockwise = hole, clockwise = normal boundary}
##' \item{\pkg{spatstat}:}{do \emph{not repeat} first vertex,
##' anticlockwise = normal boundary, clockwise = hole. This convention is also
##' used in \code{xylist}.}
##' \item{\pkg{gpclib}:}{Unfortunately, there seems to be no convention
##' for the specification of polygons of class \code{"gpc.poly"}.}
##' }
##'
##' @param object an object of one of the supported spatial classes.
##' @param ... (unused) argument of the generic.
##' @return Applying \code{xylist} to a polygon object, one gets a simple list,
##' where each component (polygon) is a list of \code{"x"} and \code{"y"}
##' coordinates. These represent vertex coordinates following \pkg{spatstat}'s
##' \code{"owin"} convention (anticlockwise order without repeating any vertex).
##' The opposite vertex order can be retained for the \pkg{sp}-classes
##' by the non-default use with \code{reverse=FALSE}.
##' @author Sebastian Meyer
##' @name xylist
##' @keywords spatial methods
##' @export
xylist <- function (object, ...) UseMethod("xylist")

##' @rdname xylist
##' @importFrom spatstat as.polygonal
##' @export
xylist.owin <- function (object, ...)
{
    as.polygonal(object)$bdry
}

##' @rdname xylist
##' @export
xylist.gpc.poly <- function (object, ...)
{
    xylist.owin(gpc2owin(object, check = FALSE))
}

##' @rdname xylist
##' @inheritParams xylist.Polygons
##' @export
xylist.SpatialPolygons <- function (object, reverse = TRUE, ...)
{
    unlist(lapply(object@polygons, xylist.Polygons, reverse=reverse, ...),
           recursive=FALSE, use.names=FALSE)
}

##' @rdname xylist
##' @param reverse logical (\code{TRUE}) indicating if the vertex order of the
##' \pkg{sp} classes should be reversed to get the \code{xylist}/\code{owin}
##' convention.
##' @import sp
##' @export
xylist.Polygons <- function (object, reverse = TRUE, ...)
{
    lapply(object@Polygons, function (sr) {
        coords <- coordinates(sr)
        n <- nrow(coords) - 1L   # number of vertices
        idxs <- if (reverse) seq.int(n,1) else seq_len(n)
        list(x = coords[idxs,1L], y = coords[idxs,2L])
             #area = sr@area, hole = sr@hole
    })
}

##' @rdname xylist
##' @import methods
##' @export
xylist.Polygon <- function (object, reverse = TRUE, ...)
    xylist.Polygons(as(object,"Polygons"), reverse=reverse, ...)

##' @rdname xylist
##' @importFrom grDevices xy.coords
##' @export
xylist.default <- function (object, ...) {
    lapply(object, function (xy) {
        poly <- xy.coords(xy)[c("x","y")]
        if (isClosed(poly)) {
            sel <- seq_len(length(poly$x) - 1L)
            poly$x <- poly$x[sel]
            poly$y <- poly$y[sel]
        }
        poly
    })
}
