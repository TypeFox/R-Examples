#--------------------------------------------------------------------
#   simage.R (npsp package)
#--------------------------------------------------------------------
#   simage  S3 generic
#       simage.default
#       simage.data.grid
#   plot.np.den
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
# simage 
#--------------------------------------------------------------------
#' Image plot with a color scale
#'
#' \code{simage} (generic function) draws an image (a grid of colored rectangles) 
#' and (optionally) adds a legend strip with the color scale 
#' (calls \code{\link{splot}} and \code{\link{image}}). 
#'
#' @seealso \code{\link{splot}}, \code{\link{spoints}}, \code{\link{spersp}}, 
#' \code{\link{image}}, \code{\link[fields]{image.plot}}, \code{\link{data.grid}}.
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
simage <- function(x, ...) UseMethod("simage")
# S3 generic function simage
#--------------------------------------------------------------------



#--------------------------------------------------------------------
# simage.default
#--------------------------------------------------------------------
#' @rdname simage  
#' @method simage default
#' @param x grid values for \code{x} coordinate. If \code{x} is a list, 
#' its components \code{x$x} and \code{x$y} are used for \code{x}  
#' and \code{y}, respectively. For compatibility with \code{\link{image}}, if the 
#' list has component \code{z} this is used for \code{s}.
#' @param y grid values for \code{y} coordinate.
#' @param s matrix containing the values to be used for coloring the rectangles (NAs are allowed). 
#' Note that \code{x} can be used instead of \code{s} for convenience.
#' @param legend logical; if \code{TRUE} (default), the plotting region is splitted into two parts,
#' drawing the image plot in one and the legend with the color scale in the other.
#' If \code{FALSE} only the image plot is drawn and the arguments related 
#' to the legend are ignored (\code{\link{splot}} is not called).
#' @param ... additional graphical parameters (to be passed to \code{\link{image}} 
#' or \code{simage.default}; e.g. \code{xlim, ylim,} ...). NOTE:
#' graphical arguments passed here will only have impact on the main plot. 
#' To change the graphical defaults for the legend use the \code{\link{par}} 
#' function beforehand (e.g. \code{par(cex.lab = 2)} to increase colorbar labels). 
#' @return Invisibly returns a list with the following 3 components:
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
#' simage( x1, x2, trend, main = 'Trend')
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
#' simage( x1, x2, y, main = 'Data')
#' simage(lp, main = 'Estimated trend')
#' par(old.par)
#' @export
#--------------------------------------------------------------------
simage.default <- function(x = seq(0, 1, len = nrow(s)), y = seq(0, 1, 
    len = ncol(s)), s, slim = range(s, finite = TRUE), col = jet.colors(128), 
    breaks = NULL, legend = TRUE, horizontal = FALSE, legend.shrink = 0.8, 
    legend.width = 1.2, legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
    bigplot = NULL, smallplot = NULL, lab.breaks = NULL, axis.args = NULL, 
    legend.args = NULL, graphics.reset = FALSE, xlab = NULL, ylab = NULL,
    ...) {
#--------------------------------------------------------------------
    if (missing(s)) {
        if (!missing(x)) {
            if (is.list(x)) {
                s <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                s <- x
                if (!is.matrix(s)) 
                    stop("argument 's' must be a matrix")                
                x <- seq.int(0, 1, length.out = nrow(s))
            }
        }
        else stop("no 's' matrix specified")
    }
    else if (is.list(x)) {
        xn <- deparse(substitute(x))
        if (missing(xlab)) xlab <- paste(xn, "x", sep = "$")
        if (missing(ylab)) ylab <- paste(xn, "y", sep = "$")
        y <- x$y
        x <- x$x
    }
    if (!is.matrix(s))
        if (missing(x) | missing(y)) stop("argument 's' must be a matrix")
        else dim(s) <- c(length(x), length(y))     
    if (is.null(xlab)) 
        xlab <- if (!missing(x)) 
            deparse(substitute(x))
        else "X"
    if (is.null(ylab)) 
        ylab <- if (!missing(y)) 
            deparse(substitute(y))
        else "Y"
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
    }
    image(x, y, s, xlab = xlab, ylab = ylab, col = col, breaks = breaks, ...)
    box()   
    if (graphics.reset) par(res$old.par)    
    return(invisible(res))        
#--------------------------------------------------------------------
}   # simage.default



#--------------------------------------------------------------------
#' @rdname simage  
#' @method simage data.grid
#' @param data.ind integer or character with the index or name of the component 
#' containing the values to be used for coloring the rectangles.
#' @export
simage.data.grid <- function(x, data.ind = 1, xlab = NULL, ylab = NULL, ...) {
#--------------------------------------------------------------------
    if (!inherits(x, "data.grid") | x$grid$nd != 2L)
        stop("function only works for two-dimensional gridded data ('data.grid'-class objects)")
    coorvs <- coordvalues(x)
    ns <- names(coorvs)
    if (is.null(xlab)) xlab <- ns[1]
    if (is.null(ylab)) ylab <- ns[2]
    res <- simage.default(coorvs[[1]], coorvs[[2]], s = x[[data.ind]],  
        xlab = xlab, ylab = ylab, ...)  
    return(invisible(res))
#--------------------------------------------------------------------
} # simage.grid.par



#--------------------------------------------------------------------
#' @rdname simage 
#' @method plot np.den 
#' @description \code{plot.np.den} calls \code{simage.data.grid} 
#' (\code{\link{contour}} and \code{\link{points}} also by default). 
#' @param log logical; if \code{TRUE} (default), \code{log(x$est)} is ploted.
#' @param contour logical; if \code{TRUE} (default), contour lines are added.
#' @param points logical; if \code{TRUE} (default), points at \code{x$data$x} are drawn.
#' @export
plot.np.den <- function(x, y = NULL, log = TRUE, contour = TRUE, points = TRUE, 
                    col = hot.colors(128), ...){
#    if (!inherits(x, "data.grid") | x$grid$nd != 2L)
#        stop("function only works for two-dimensional gridded data ('data.grid'-class objects)")
    if (log) x$est <- log(x$est)
    ret <- simage(x, col = col, ...)
    if (contour) contour(x, add = TRUE)
    if (points) points(x$data$x, pch = 21, bg = 'black', col = 'darkgray' )
    return(invisible(ret))
#--------------------------------------------------------------------
} # plot.np.den