#--------------------------------------------------------------------
#   npsp-plot.R (npsp package)
#   S3 and S4 plot methods for npsp-objects
#--------------------------------------------------------------------
#   svar.plot S3 plot methods
#     plot.fitsvar(x, y, legend, xlab, ylab, ylim, lwd, add, ...)
#     plot.svar.bin(x, y, xlab, ylab, ylim, add, ...)
#     plot.np.svar(x, y, xlab, ylab, ylim, add, ...)
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2014
#--------------------------------------------------------------------
# PENDENTE:
#   - @examples
#   - curve(sv(svm, x), add = TRUE)
#--------------------------------------------------------------------


#' @name svar.plot
#' @title Plot a semivariogram object
#' @description Utilities for plotting pilot semivariograms or fitted models.
#' @seealso
#' \code{\link{svariso}}, \code{\link{np.svariso}}, \code{\link{fitsvar.sb.iso}}.
NULL


#--------------------------------------------------------------------
#' @rdname svar.plot
#' @method plot fitsvar
#' @description \code{plot.fitsvar} plots a fitted variogram model.
#' @param x a variogram object. Typically an output of functions
#'   \code{\link{np.svariso}} or \code{\link{fitsvar.sb.iso}}.
#' @param y ignored argument.
#' @param legend logical; if \code{TRUE} (default), a legend is added to the plot.
#' @param xlab label for the x axis (defaults to "distance").
#' @param ylab label for the y axis (defaults to "semivariance").
#' @param ylim y-limits.
#' @param lwd line widths for points (estimates) and lines (fitted model) respectively.
#' @param ...  additional graphical parameters (see \code{\link{par}}).
#' @export
plot.fitsvar <- function(x, y = NULL, legend = TRUE, xlab = "distance", ylab = "semivariance",
                          ylim = c(0, 1.25*max(x$fit$sv, na.rm = TRUE)), lwd = c(1, 2), add = FALSE, ...) {
    if (!inherits(x, "isotropic"))
      stop("function only works for isotropic variograms.")
    if (add)   
        with(x$fit, lines(u, fitted.sv, ... ))
    else {
        with(x$fit, plot(u, sv, xlab = xlab, ylab = ylab, ylim = ylim, lwd = lwd[1], ...))
        with(x$fit, lines(u, fitted.sv, lwd = lwd[2]))
        if (legend) legend("bottomright", legend = c("estimates", "fitted model"),
            lty = c(NA, 1), pch = c(1, NA), lwd = lwd)
    } 
    invisible(NULL)       
#--------------------------------------------------------------------
} # plot.fitsvar


#--------------------------------------------------------------------
#' @rdname svar.plot
#' @method plot svar.bin
#' @description \code{plot.svar.bin} plots the binned semivariances.
#' @param add logical; if \code{TRUE} the semivariogram plot is just added 
#' to the existing plot.
#' @export
plot.svar.bin <- function(x, y = NULL, xlab = "distance", ylab = "semivariance",
                          ylim = c(0,max(x$biny, na.rm = TRUE)), add = FALSE, ...) {
    if (x$grid$nd != 1L)
        stop("function only works for isotropic variograms.")
    if (add)   
        points(coords(x), x$biny, ... )
    else
        plot(coords(x), x$biny, xlab = xlab, ylab = ylab, ylim = ylim, ... )
#--------------------------------------------------------------------
} # plot.svar.bin
        


#--------------------------------------------------------------------
#' @rdname svar.plot
#' @method plot np.svar
#' @description \code{plot.np.svar} plots a local polynomial estimate of the semivariogram.
#' @export
plot.np.svar <- function(x, y = NULL, xlab = "distance", ylab = "semivariance",
                          ylim = c(0,max(x$biny, na.rm = TRUE)), add = FALSE, ...) {
    if (x$grid$nd != 1L)
        stop("function only works for isotropic variograms.")
    if (add)   
        lines(coords(x), x$est, ... )
    else {
        plot(coords(x), x$biny, xlab = xlab, ylab = ylab, ylim = ylim, ... )
        lines(coords(x), x$est, lwd = 2 )
    }    
#--------------------------------------------------------------------
} # plot.svar.bin
