################################################################################
### Part of the R package "polyCub".
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Copyright (C) 2009-2015 Sebastian Meyer
### Time-stamp: <[polyCub.midpoint.R] 2015-02-25 20:53 (CET) by SM>
################################################################################


#' Two-Dimensional Midpoint Rule
#'
#' The surface is converted to a binary pixel image
#' using the \code{\link[spatstat]{as.im.function}} method from package
#' \pkg{spatstat} (Baddeley and Turner, 2005).
#' The integral under the surface is then approximated as the
#' sum over (pixel area * f(pixel midpoint)).
#' 
#' @inheritParams plotpolyf
#' @param polyregion a polygonal integration domain.
#' It can be any object coercible to the \pkg{spatstat} class
#' \code{"\link[spatstat]{owin}"} via a corresponding
#' \code{\link[spatstat]{as.owin}}-method.
#' Note that this includes polygons of the classes \code{"gpc.poly"} and
#' \code{"\linkS4class{SpatialPolygons}"}, because \pkg{polyCub} defines
#' methods \code{\link{as.owin.gpc.poly}} and
#' \code{\link{as.owin.SpatialPolygons}}, respectively.
#' @param eps width and height of the pixels (squares),
#' see \code{\link[spatstat]{as.mask}}.
#' @param dimyx number of subdivisions in each dimension,
#' see \code{\link[spatstat]{as.mask}}.
#' @param plot logical indicating if an illustrative plot of the numerical
#' integration should be produced.
#' @return The approximated value of the integral of \code{f} over
#' \code{polyregion}.
#' @references
#' Baddeley, A. and Turner, R. (2005).
#' \pkg{spatstat}: an \R package for analyzing spatial point patterns.
#' \emph{Journal of Statistical Software}, \bold{12} (6), 1-42.
#' @keywords math spatial
#' @family polyCub-methods
#' @import sp
#' @importFrom spatstat as.im.function plot.im integral.im
#' @importFrom grDevices gray
#' @examples # see example(polyCub)
#' @export
## NOTE: we don't import graphics::plot since it is already imported via sp

polyCub.midpoint <- function (polyregion, f, ...,
                              eps = NULL, dimyx = NULL, plot = FALSE)
{
    ## as.im needs seperate x and y arguments
    fxy <- function (x, y, ...) f(cbind(x,y), ...)

    ## calculate pixel values of fxy
    IM <- tryCatch(
          as.im.function(X=fxy, W=polyregion, ..., eps=eps, dimyx=dimyx),
          error = function (e) {
              ## if eps was to small such that the dimensions of the image would
              ## be too big then the operation matrix(TRUE, nr, nc) throws an
              ## error. (try e.g. devnull <- matrix(TRUE, 1e6,1e6))
              ## unfortunately, it is not clear what we should do in this
              ## case... => stop
              stop("inapplicable choice of bandwidth (eps=", format(eps),
                   ") in midpoint rule:\n", e)
          })
    
### ILLUSTRATION ###
    if (plot) {
        plot.im(IM, axes=TRUE, col=gray(31:4/35), main="")
        ## add evaluation points
        #with(IM, points(expand.grid(xcol, yrow), col=!is.na(v), cex=0.5))
        plot(polyregion, add=TRUE, poly.args=list(lwd=2), lwd=2)
        ##<- two 'lwd'-specifications such that it works with owin and gpc.poly
    }
####################
    
    ## return the approximated integral
    integral.im(IM)
}
