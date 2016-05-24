#--------------------------------------------------------------------
#   spersp.R (npsp package)
#--------------------------------------------------------------------
#   spersp  S3 generic
#       spersp.default
#       spersp.data.grid
#       persp.data.grid
#
#   Based on image.plot and drape.plot functions from package fields:
#   fields, Tools for spatial data
#   Copyright 2004-2013, Institute for Mathematics Applied Geosciences
#   University Corporation for Atmospheric Research
#   Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
#   (c) R. Fernandez-Casal         Last revision: Mar 2014
#--------------------------------------------------------------------



#--------------------------------------------------------------------
# spersp 
#--------------------------------------------------------------------
#' Perspective plot with a color scale
#'
#' \code{spersp} (generic function) draws a perspective plot of a surface over 
#' the \code{x-y} plane with the facets being filled with different colors 
#' and (optionally) adds a legend strip with the color scale 
#' (calls \code{\link{splot}} and \code{\link{persp}}). 
#'
#' @seealso \code{\link{splot}}, \code{\link{spoints}}, \code{\link{simage}}, 
#' \code{\link{image}}, \code{\link[fields]{image.plot}}, \code{\link{data.grid}}, 
#' \code{\link{persp}}.
#' @section Side Effects: After exiting, the plotting region may be changed 
#' (\code{\link{par}("plt")}) to make it possible to add more features to the plot
#' (set \code{graphics.reset = FALSE} to avoid this).
#' @author
#' Based on \code{\link[fields]{image.plot}} function from package \pkg{fields}:
#' fields, Tools for spatial data. 
#' Copyright 2004-2013, Institute for Mathematics Applied Geosciences. 
#' University Corporation for Atmospheric Research.
#' 
#' Modified by Ruben Fernandez-Casal <rubenfcasal@@gmail.com>. 
#' @keywords hplot
#' @export
#--------------------------------------------------------------------
spersp <- function(x, ...) UseMethod("spersp")
# S3 generic function spersp
#--------------------------------------------------------------------



#--------------------------------------------------------------------
# spersp.default
#--------------------------------------------------------------------
#' @rdname spersp  
#' @method spersp default
#' @param x grid values for \code{x} coordinate. If \code{x} is a list, 
#' its components \code{x$x} and \code{x$y} are used for \code{x}  
#' and \code{y}, respectively. If the list has component \code{z} this is used 
#' for \code{z}.
#' @param y grid values for \code{y} coordinate.
#' @param z matrix containing the values to be plotted (NAs are allowed).  
#' Note that \code{x} can be used instead of \code{z} for convenience.
#' @param s matrix containing the values used for coloring the facets. 
#' @param legend logical; if \code{TRUE} (default), the plotting region is splitted into two parts,
#' drawing the perspective plot in one and the legend with the color scale in the other.
#' If \code{FALSE} only the (coloured) perspective plot is drawn and the arguments related 
#' to the legend are ignored (\code{\link{splot}} is not called).
#' @param zlab label for the z axis, defaults to a description of \code{z}.
#' @param theta x-y rotation angle for perspective (azimuthal direction).
#' @param phi z-angle for perspective (colatitude). 	
#' @param ticktype character; \code{"simple"} draws just an arrow parallel to the axis 
#' to indicate direction of increase; \code{"detailed"} draws normal ticks as per 2D plots.
#' @param cex.axis magnification to be used for axis annotation (relative to the 
#' current setting of \code{\link{par}("cex")}).
#' @param ... additional graphical parameters (to be passed to \code{\link{persp}} 
#' or \code{spersp.default}; e.g. \code{xlim, ylim, zlim,} ...). NOTE:
#' graphical arguments passed here will only have impact on the main plot. 
#' To change the graphical defaults for the legend use the \code{\link{par}} 
#' function beforehand (e.g. \code{par(cex.lab = 2)} to increase colorbar labels). 
#' @return Invisibly returns a list with the following 4 components:
#' \item{pm}{the viewing transformation matrix (see \code{\link{persp}} for details), 
#' a 4 x 4 matrix that can be used to superimpose additional graphical elements 
#' using the function \code{\link{trans3d}}.} 
#' \item{bigplot}{plot coordinates of the main plot. These values may be useful for 
#' drawing a plot without the legend that is the same size as the plots with legends.}
#' \item{smallplot}{plot coordinates of the secondary plot (legend strip).}
#' \item{old.par}{previous graphical parameters (\code{par(old.par)} 
#' will reset plot parameters to the values before entering the function).}
#' @inheritParams splot
#' @inheritParams spoints
#' @examples
#' 
#' #
#' # Regularly spaced 2D data
#' nx <- c(40, 40) # ndata =  prod(nx)
#' x1 <- seq(-1, 1, length.out = nx[1])
#' x2 <- seq(-1, 1, length.out = nx[2])
#' trend <- outer(x1, x2, function(x,y) x^2 - y^2) 
#' spersp( x1, x2, trend, main = 'Trend', zlab = 'y')
#' 
#' #
#' # Multiple plots 
#' set.seed(1)
#' y <- trend + rnorm(prod(nx), 0, 0.1)
#' x <- as.matrix(expand.grid(x1 = x1, x2 = x2)) # two-dimensional grid
#' # local polynomial kernel regression
#' lp <- locpol(x, y, nbin = nx, h =  diag(c(0.3, 0.3)))
#' # 1x2 plot
#' old.par <- par(mfrow = c(1,2))
#' spersp( x1, x2, y, main = 'Data')
#' spersp(lp, main = 'Estimated trend', zlab = 'y')
#' par(old.par)
#' @export
#--------------------------------------------------------------------
spersp.default <- function(x = seq(0, 1, len = nrow(z)), y = seq(0, 1, 
    len = ncol(z)), z , s = z, slim = range(s, finite = TRUE), col = jet.colors(128), 
    breaks = NULL, legend = TRUE, horizontal = FALSE, legend.shrink = 0.8, 
    legend.width = 1.2, legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
    bigplot = NULL, smallplot = NULL, lab.breaks = NULL, axis.args = NULL, 
    legend.args = NULL, graphics.reset = FALSE, xlab = NULL, ylab = NULL,
    zlab = NULL, theta = 40, phi = 20, ticktype = "detailed", cex.axis = 0.75, ...) {
# if s is passed ( values for coloring facets ) use it
# if not use the z matrix that is also used to draw the
# perspective plot.
#--------------------------------------------------------------------
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                if (!is.matrix(z))
                    stop("argument 'z' must be a matrix")
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    } else if (is.list(x)) {
        xn <- deparse(substitute(x))
        if (missing(xlab)) xlab <- paste(xn, "x", sep = "$")
        if (missing(ylab)) ylab <- paste(xn, "y", sep = "$")
        y <- x$y
        x <- x$x
    }
    if (is.null(xlab)) 
        xlab <- if (!missing(x)) 
            deparse(substitute(x))
        else "X"
    if (is.null(ylab)) 
        ylab <- if (!missing(y)) 
            deparse(substitute(y))
        else "Y"
    if (is.null(zlab)) 
        zlab <- if (!missing(z)) 
            deparse(substitute(z))
        else "Z"
    if (!is.matrix(z))
        if (missing(x) | missing(y)) stop("argument 'z' must be a matrix")
        else dim(z) <- c(length(x), length(y)) 

    if (!missing(s) & !identical(dim(z), dim(s)))
              stop("'s' matrix dimensions must match 'z'")
    # do average before? (now analogous to 'drape.plot'...)
    if (legend) 
        # image in splot checks breaks and other parameters...
        res <- splot(slim = slim, col = col, breaks = breaks, horizontal = horizontal, 
            legend.shrink = legend.shrink, legend.width = legend.width,
            legend.mar = legend.mar, legend.lab = legend.lab,
            bigplot = bigplot, smallplot = smallplot, lab.breaks = lab.breaks,        
            axis.args = axis.args, legend.args = legend.args)
    else {
        old.par <- par(no.readonly = TRUE)
        # par(xpd = FALSE)
        res <- list(bigplot = old.par$plt, smallplot = NA, old.par = old.par)    
    }        
    if (is.null(breaks)) {
        # Compute breaks (in 'cut.default' style...)
        ds <- diff(slim)
        if (ds == 0) ds <- abs(slim[1L])
        breaks <- seq.int(slim[1L] - ds/1000, slim[2L] + ds/1000, length.out = length(col) + 1)
        # Only if !missing(slim) else breaks <- length(col) + 1?
    }
    # average s value for a facet
    nx <- nrow(s)
    ny <- ncol(s)          
    s <- 0.25 * (s[-nx,-ny] + s[-1,-ny] + s[-nx,-1] + s[-1,-1])
    # set colors
    icol <- cut(as.numeric(s), breaks, labels = FALSE, include.lowest = TRUE, right = FALSE) # Use .bincode instead of cut?      
    # call persp
    pm <- persp(x, y, z,  xlab = xlab, ylab = ylab, zlab = zlab, theta = theta, 
           phi = phi, col = col[icol], ticktype = ticktype, cex.axis = cex.axis, ...)   
    if (graphics.reset) par(res$old.par)    
    return(invisible( c(list(pm = pm), res) ))        
#--------------------------------------------------------------------
}   # spersp.default



#--------------------------------------------------------------------
#' @rdname spersp  
#' @method spersp data.grid
#' @param data.ind integer or character with the index or name of the component 
#' containing the \code{z} values to be plotted.
#' @export
spersp.data.grid <- function(x, data.ind = 1, s = x[[data.ind]], 
        xlab = NULL, ylab = NULL, zlab = NULL, ...) {
#--------------------------------------------------------------------
    if (!inherits(x, "data.grid") | x$grid$nd != 2L)
        stop("function only works for two-dimensional gridded data ('data.grid'-class objects)")
    coorvs <- coordvalues(x)
    ns <- names(coorvs)
    if (is.null(xlab)) xlab <- ns[1]
    if (is.null(ylab)) ylab <- ns[2]
    if (is.null(zlab)) 
        zlab <- if (is.character(data.ind)) data.ind else names(x)[data.ind]
    res <- spersp.default(coorvs[[1]], coorvs[[2]], z = x[[data.ind]], s = s, 
        xlab = xlab, ylab = ylab, zlab = zlab, ...)  
    return(invisible(res))
#--------------------------------------------------------------------
} # spersp.grid.par

