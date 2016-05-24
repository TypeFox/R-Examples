#--------------------------------------------------------------------
#   rgraphics.R (npsp package)
#--------------------------------------------------------------------
#   data.grid S3 methods
#     image.data.grid(x, data.ind, xlab, ylab, ...)
#     persp.data.grid(x, data.ind, xlab, ylab, zlab, ...) 
#     contour.data.grid(x, data.ind, filled, xlab, ylab, ...)
#
#   (c) R. Fernandez-Casal         Last revision: Mar 2014
#--------------------------------------------------------------------
# PENDENTE:
#   - plot.data.grid
#--------------------------------------------------------------------


#' @name rgraphics
#' @title R Graphics for gridded data
#' @description Draw an image, perspective, contour or filled contour plot for data
#' on a bidimensional regular grid (S3 methods for class "code{\link{data.grid}}").
#' @seealso \code{\link{image}}, \code{\link{persp}}, \code{\link{contour}},
#' \code{\link{filled.contour}}, \code{\link{data.grid}}.
#' @examples
#' # Regularly spaced 2D data
#' grid <- grid.par(n = c(50, 50), min = c(-1, -1), max = c(1, 1))
#'
#' f2d <- function(x) x[1]^2 - x[2]^2
#' trend <- apply(coords(grid), 1, f2d)
#' set.seed(1)
#' y <- trend + rnorm(prod(dim(grid)), 0, 0.1)
#' gdata <- data.grid(trend = trend, y = y, grid = grid)
#'
#' # perspective plot
#' persp(gdata, main = 'Trend', theta = 40, phi = 20, ticktype = "detailed")
#'
#' # filled contour plot
#' contour(gdata, main = 'Trend', filled = TRUE, color.palette = jet.colors)
#'
#' # Multiple plots with a common legend:
#' scale.range <- c(-1.2, 1.2)
#' scale.color <- jet.colors(64)
#' # 1x2 plot with some room for the legend...
#' old.par <- par(mfrow = c(1,2), omd = c(0.05, 0.85, 0.05, 0.95))
#' image(gdata, zlim = scale.range, main = 'Trend', col = scale.color)
#' contour(gdata, add = TRUE)
#' image(gdata, 'y', zlim = scale.range, main = 'Data', col = scale.color)
#' contour(gdata, 'y', add = TRUE)
#' par(old.par)
#' # the legend can be added to any plot...
#' splot(slim = scale.range, col = scale.color, add = TRUE)
NULL


#--------------------------------------------------------------------
#' @rdname rgraphics
#' @method image data.grid
#' @param x a "code{\link{data.grid}}"-class object.
#' @param data.ind integer or character with the index or name of the component
#' containing the values to be used for coloring the rectangles.
#' @param xlab label for the x axis, defaults to \code{dimnames(x)[1]}.
#' @param ylab label for the y axis, defaults to \code{dimnames(x)[2]}.
#' @param zlab label for the z axis, defaults to \code{names(x)[data.ind]}.
#' @param ... additional graphical parameters (to be passed to main plot function).
#' @export
image.data.grid <- function(x, data.ind = 1, xlab = NULL, ylab = NULL, ...) {
#--------------------------------------------------------------------
    if (!inherits(x, "data.grid") | x$grid$nd != 2L)
        stop("function only works for two-dimensional gridded data ('data.grid'-class objects)")
    coorvs <- coordvalues(x)
    ns <- names(coorvs)
    if (is.null(xlab)) xlab <- ns[1]
    if (is.null(ylab)) ylab <- ns[2]
    image(coorvs[[1]], coorvs[[2]], z = x[[data.ind]],
        xlab = xlab, ylab = ylab, ...)
    return(invisible())
#--------------------------------------------------------------------
} # simage.grid.par



#--------------------------------------------------------------------
#' @rdname rgraphics
#' @method persp data.grid
#' @export
persp.data.grid <- function(x, data.ind = 1,
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
    res <- persp(coorvs[[1]], coorvs[[2]], z = x[[data.ind]],
        xlab = xlab, ylab = ylab, zlab = zlab, ...)
    return(invisible(res))
#--------------------------------------------------------------------
} # spersp.grid.par



#--------------------------------------------------------------------
#' @rdname rgraphics
#' @method contour data.grid
#' @param filled logical; if \code{FALSE} (default), function \code{\link{contour}}
#' is called, otherwise \code{\link{filled.contour}}.
#' @export
contour.data.grid <- function(x, data.ind = 1, filled = FALSE, xlab = NULL, ylab = NULL, ...) {
#--------------------------------------------------------------------
    if (!inherits(x, "data.grid") | x$grid$nd != 2L)
        stop("function only works for two-dimensional gridded data ('data.grid'-class objects)")
    coorvs <- coordvalues(x)
    ns <- names(coorvs)
    if (is.null(xlab)) xlab <- ns[1]
    if (is.null(ylab)) ylab <- ns[2]
    if (filled)
        filled.contour(coorvs[[1]], coorvs[[2]], z = x[[data.ind]],
            xlab = xlab, ylab = ylab, ...)
    else
        contour(coorvs[[1]], coorvs[[2]], z = x[[data.ind]],
            xlab = xlab, ylab = ylab, ...)
    return(invisible())
#--------------------------------------------------------------------
} # contour.grid.par


