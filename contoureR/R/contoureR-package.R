#' contoureR: Contouring of Non-Regular Three-Dimensional Data
#' 
#' Create contour lines for a non regular series of points, potentially from a non-regular canvas.
#' 
#' The \code{contoureR} package executes linear interpolation on a delaunay triangulated 
#' mesh strung between three-dimensional (3D) points supplied by the user. Contours are calculated 
#' across the surface constrained by the convex hull of the supplied data. 
#' 
#' Usually, the well known functions such as \code{\link{contourLines}} from the \code{\link{grDevices}} 
#' package, expect (or rather, require) data to be regular, this means that a rectangular array or matrix of 
#' \code{x} and \code{y} coordinate pairs, each with a corresponding \code{z} value is to be modelled -- that
#' is to say the cartesian product of a numeric vector of \code{x} values of length \code{n}, 
#' with a numeric vector of \code{y} values having length \code{m}, used to produce a set of
#' \code{(m x n)} unique points that have been concurrently provided with \code{exactly (m x n) z values}. 
#' 
#' By restricting values to the above format, this in turn limits the region of analysis to square/rectangular 
#' canvasses (ie plane defined by geometric and orthogonal vectors parallel to the \code{x} and \code{y} axes and range bound by the \code{[xmin,xmax]} 
#' and \code{[ymin,ymax]} in the above \code{x} and \code{y} input numeric vectors, respectively). 
#' This restriction, from time-to-time, can be very inconvenient, and is a primary objective and purpose for
#' the creation of this package.
#' 
#' As suggested in the previous paragraph, the \code{contoureR} package, on the other hand, has no such 
#' orthogonality / regularity requirement and can therefore be applied over obscurely shaped regions such as 
#' triangles, circles, polygons and the like. To demonstrate this, in the example 
#' provided on the current page, an equation is contoured, where firstly the \code{x} and \code{y} data is 
#' randomly selected (non regular), and then the set of values is subsequently constrained by a 
#' bounding (limiting) circle. 
#' 
#' Note, for the moment, the only restriction is that for polygon-type regions to be
#' modelled, then these regions must \strong{not} have holes, since these will be filled coarsely when the 
#' Deleaunaymesh gets generated, however, in future revisions, this obstacle should be easily addressed via 
#' parameter defining a manual exclusion list of points.
#' @examples
#' # Contour Lines for a Function, Constrained to a limited domain
#' # Example of the provision of non-regular data
#' library(contoureR)
#' library(ggplot2)
#' a  = -2; b = +2; n  = 150
#' x  = runif(n*n,a,b)
#' y  = runif(n*n,a,b)
#' df = data.frame(x,y)
#' df$z   = with(df,-x*y*exp(-x^2-y^2))
#' df.sub = subset(df,x^2 + y^2 < 2)
#' df.cnt = getContourLines(df.sub,nlevels=20)
#' ggplot(data=df.cnt,aes(x,y,group=Group,colour=z)) + geom_path() + theme_bw()
#' @seealso \code{\link{getContourLines}}, \code{\link{contourLinesR}} and \code{\link{contourWalker}}
#' @rdname contoureR
#' @name contoureR
NULL